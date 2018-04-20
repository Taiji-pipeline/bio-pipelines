{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE OverloadedLists      #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE ScopedTypeVariables  #-}
{-# LANGUAGE TemplateHaskell      #-}
{-# LANGUAGE TypeFamilies         #-}
{-# LANGUAGE TypeOperators        #-}
{-# LANGUAGE UndecidableInstances #-}
module Bio.Pipeline.NGS.Utils where

import           Bio.Data.Bam                (bamToBed, readBam,
                                              sortedBamToBedPE, withBamFile)
import           Bio.Data.Bed                (BED, BED3, BEDConvert (..),
                                              BEDLike (..))
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib           (gzip, ungzip)
import           Data.Maybe                  (fromJust)
import           Data.Promotion.Prelude.List (Elem)
import           Data.Singletons             (SingI)
import qualified Data.Text                   as T
import qualified Data.Text.IO                as T
import           Shelly                      (escaping, fromText, mv, run_,
                                              shelly, silently)
import           System.IO.Temp              (withTempDirectory)

-- | Remove low quality and redundant tags, fill in mate information.
filterBam_ :: (SingI tags, tags' ~ (Insert' 'Sorted tags))
           => FilePath  -- ^ output
           -> File tags 'Bam
           -> IO (File tags' 'Bam)
filterBam_ output fl = withTempDirectory "./" "tmp_filt_dir." $ \tmp -> do
    let input = T.pack $ fl^.location
    shelly $ escaping False $ silently $ do
        let tmp_filt = T.pack $ tmp ++ "/tmp_filt.bam"
            tmp_fixmate = T.pack $ tmp ++ "/tmp_fixmate.bam"
            tmp_sort = T.pack $ tmp ++ "/tmp_sort"
        run_ "samtools" $ ["view"] ++
            (if isPair then ["-f", "2"] else []) ++
            ["-F", "0x70c", "-q", "30", "-u", input] ++
            ( if isPair
                then [ "|", "samtools", "sort", "-", "-n", "-T", tmp_sort
                    , "-l", "0", "-o", tmp_filt ]
                else [ "|", "samtools", "sort", "-", "-T", tmp_sort
                    , "-l", "9", "-o", T.pack output ] )
        when isPair $ do
            run_ "samtools" ["fixmate", "-r", tmp_filt, tmp_fixmate]
            run_ "samtools" [ "view", "-F", "1804", "-f", "2", "-u"
                , tmp_fixmate, "|", "samtools", "sort", "-", "-T"
                , tmp_sort, "-l", "9", "-o", T.pack output ]

    return $ location .~ output $ emptyFile
  where
    isPair = fl `hasTag` PairedEnd

-- | Remove duplicates
removeDuplicates_ :: SingI tags
                  => FilePath
                  -> FilePath
                  -> File tags 'Bam
                  -> IO (File tags 'Bam)
removeDuplicates_ picardPath output input =
    withTempDirectory "./" "tmp_picard_dir." $ \tmp -> shelly $ do
        let qcFile = tmp ++ "/picard.qc"
            markdupTmp = tmp++"/dup_marked.bam"
            filtTmp = tmp++"/dup_filt.bam"
        -- Mark duplicates
        run_ "java" ["-Xmx4G", "-jar", T.pack picardPath
            , "MarkDuplicates", T.pack $ "INPUT=" ++ (input^.location)
            , T.pack $ "OUTPUT=" ++ markdupTmp
            , T.pack $ "TMP_DIR=" ++ tmp
            , T.pack $ "METRICS_FILE=" ++ qcFile
            , "VALIDATION_STRINGENCY=LENIENT"
            , "ASSUME_SORT_ORDER=coordinate", "REMOVE_DUPLICATES=false"]

        -- Remove duplicates.
        escaping False $ run_ "samtools" [ "view", "-F", "0x70c", "-b"
            , T.pack markdupTmp, ">", T.pack filtTmp ]

        -- Re-sort by names for pairedend sequencing
        if input `hasTag` PairedEnd
            then run_ "samtools" [ "sort", T.pack filtTmp, "-n", "-T"
                , T.pack $ tmp ++ "/tmp_sort", "-o", T.pack output ]
            else mv (fromText $ T.pack filtTmp) $ fromText $ T.pack output

        qc <- liftIO $ T.readFile qcFile
        return $ info .~ [("QC", qc)] $ location .~ output $ emptyFile

bam2Bed_ :: FilePath
         -> (BED -> Bool)  -- ^ Filtering function
         -> File tags 'Bam -> IO (File (Insert' 'Gzip tags) 'Bed)
bam2Bed_ output fn fl = do
    withBamFile (fl^.location) $ \h -> runConduit $ readBam h .| bamToBed .|
        filterC fn .| mapC toLine .| unlinesAsciiC .| gzip .| sinkFileBS output
    return $ location .~ output $ emptyFile
{-# INLINE bam2Bed_ #-}

-- | Convert name sorted BAM to BEDPE suitable for MACS2.
bam2BedPE_ :: Elem 'Sorted tags ~ 'True
           => String
           -> ((BED, BED) -> Bool)
           -> File tags 'Bam
           -> IO (File (Insert' 'Gzip tags) 'Bed)
bam2BedPE_ output fn fl = do
    withBamFile (fl^.location) $ \h -> runConduit $ readBam h .|
        sortedBamToBedPE .| filterC fn .| concatMapC f .| mapC toLine .|
        unlinesAsciiC .| gzip .| sinkFileBS output
    return $ location .~ output $ emptyFile
  where
    f (b1, b2)
        | b1^.chrom /= b2^.chrom || b1^.strand == b2^.strand = Nothing
        | otherwise =
            let left = if fromJust (b1^.strand)
                    then b1^.chromStart else b2^.chromStart
                right = if not (fromJust $ b2^.strand)
                    then b2^.chromEnd else b1^.chromEnd
            in if left < right
                  then Just (asBed (b1^.chrom) left right :: BED3)
                  else error "Left coordinate is larger than right coordinate."
{-# INLINE bam2BedPE_ #-}

-- | Merge multiple BED files.
concatBed_ :: (Elem 'Gzip tags1 ~ 'False, Elem 'Gzip tags2 ~ 'True)
           => FilePath
           -> [Either (File tags1 'Bed) (File tags2 'Bed)]
           -> IO (File '[Gzip] 'Bed)
concatBed_ output fls = do
    runResourceT $ runConduit $ source .| gzip .| sinkFile output
    return $ location .~ output $ emptyFile
  where
    source = forM_ fls $ \fl -> case fl of
        Left fl'  -> sourceFileBS (fl'^.location)
        Right fl' -> sourceFileBS (fl'^.location) .| ungzip
{-# INLINE concatBed_ #-}
