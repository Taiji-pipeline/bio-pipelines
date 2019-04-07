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

import           Bio.Data.Bam                (bamToBedC, streamBam, getBamHeader,
                                              sortedBamToBedPE )
import           Bio.Data.Bed                (BED, BED3, BEDConvert (..),
                                              BEDLike (..))
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib           (gzip, ungzip)
import           Data.Maybe                  (fromJust)
import           Data.Singletons.Prelude      (Elem, If, SingI)
import           Data.Singletons.Prelude.List (Delete)
import qualified Data.Text                   as T
import           Shelly                      (lastStderr, escaping, run, run_, shelly, silently, bash_, bashPipeFail)
import           System.IO.Temp              (withTempDirectory)

-- | Remove low quality and redundant tags, fill in mate information.
filterBam :: ( SingI tags, tags' ~ If (Elem PairedEnd tags)
               (Insert' 'NameSorted (Delete 'CoordinateSorted tags)) tags )
          => FilePath  -- ^ temp dir
          -> FilePath  -- ^ output
          -> File tags 'Bam
          -> IO (File tags' 'Bam)
filterBam tmpDir output fl = withTempDirectory tmpDir "tmp_filt_dir." $ \tmp -> do
    let input = T.pack $ fl^.location
        tmp_filt = T.pack $ tmp ++ "/tmp_filt.bam"
        tmp_fixmate = T.pack $ tmp ++ "/tmp_fixmate.bam"
        tmp_sort = T.pack $ tmp ++ "/tmp_sort"
    shelly $ escaping False $ silently $ if isPair
        then do
            bashPipeFail bash_ "samtools"
                [ "view", "-f", "2", "-F", "0x70c", "-q", "30"
                , "-u", input, "|", "samtools", "sort", "-", "-n", "-T", tmp_sort
                , "-l", "0", "-o", tmp_filt ]
            run_ "samtools" ["fixmate", "-r", "-m", tmp_filt, tmp_fixmate]
            run_ "samtools" [ "view", "-F", "1804", "-f", "2", "-u"
                , tmp_fixmate, ">", T.pack output ]
        else run_ "samtools" [ "view", "-F", "0x70c", "-q", "30", "-u", input
            , ">", T.pack output ]
    return $ location .~ output $ emptyFile
  where
    isPair = fl `hasTag` PairedEnd
{-# INLINE filterBam #-}

sortBam :: ( SingI tags, tags' ~ Insert' 'CoordinateSorted
             (Delete 'NameSorted tags) )
        => FilePath    -- ^ temp dir
        -> FilePath    -- ^ output
        -> File tags 'Bam
        -> IO (File tags' 'Bam)
sortBam tmpDir output fl = withTempDirectory tmpDir "tmp_sort_dir." $ \tmp -> do
    let input = T.pack $ fl^.location
        tmp_sort = T.pack $ tmp ++ "/tmp_sort"
    shelly $ silently $ run_ "samtools"
        [ "sort", input, "-T", tmp_sort, "-l", "9", "-o", T.pack output ]
    return $ location .~ output $ emptyFile
{-# INLINE sortBam #-}

sortBamByName :: ( SingI tags, tags' ~ Insert' 'NameSorted
                   (Delete 'CoordinateSorted tags) )
              => FilePath    -- ^ temp dir
              -> FilePath    -- ^ output
              -> File tags 'Bam
              -> IO (File tags' 'Bam)
sortBamByName tmpDir output fl = withTempDirectory tmpDir "tmp_sort_dir." $ \tmp -> do
    let input = T.pack $ fl^.location
        tmp_sort = T.pack $ tmp ++ "/tmp_sort"
    shelly $ silently $ run_ "samtools"
        [ "sort", input, "-n", "-T", tmp_sort, "-l", "9", "-o", T.pack output ]
    return $ location .~ output $ emptyFile
{-# INLINE sortBamByName #-}

-- | Remove duplicates
removeDuplicates :: Elem 'CoordinateSorted tags ~ 'True
                 => FilePath
                 -> File tags 'Bam
                 -> IO (File tags 'Bam)
removeDuplicates output input = shelly $ do
    _ <- run "samtools" [ "markdup", T.pack $ input^.location
        , T.pack output, "-r", "-s" ]
    qc <- lastStderr
    return $ info .~ [("QC", qc)] $ location .~ output $ emptyFile
{-# INLINE removeDuplicates #-}

bam2Bed :: FilePath
        -> (BED -> Bool)  -- ^ Filtering function
        -> File tags 'Bam -> IO (File (Insert' 'Gzip tags) 'Bed)
bam2Bed output fn fl = do
    header <- getBamHeader $ fl^.location
    runResourceT $ runConduit $ streamBam (fl^.location) .| bamToBedC header .|
        filterC fn .| mapC toLine .| unlinesAsciiC .| gzip .| sinkFileBS output
    return $ location .~ output $ emptyFile
{-# INLINE bam2Bed #-}

-- | Convert name sorted BAM to BEDPE suitable for MACS2.
bam2BedPE :: Elem 'NameSorted tags ~ 'True
          => String
          -> ((BED, BED) -> Bool)
          -> File tags 'Bam
          -> IO (File (Insert' 'Gzip tags) 'Bed)
bam2BedPE output fn fl = do
    header <- getBamHeader $ fl^.location
    runResourceT $ runConduit $ streamBam (fl^.location) .|
        sortedBamToBedPE header .| filterC fn .| concatMapC f .| mapC toLine .|
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
{-# INLINE bam2BedPE #-}

-- | Merge multiple BED files.
concatBed :: (Elem 'Gzip tags1 ~ 'False, Elem 'Gzip tags2 ~ 'True)
          => FilePath
          -> [Either (File tags1 'Bed) (File tags2 'Bed)]
          -> IO (File '[Gzip] 'Bed)
concatBed output fls = do
    runResourceT $ runConduit $ source .| gzip .| sinkFile output
    return $ location .~ output $ emptyFile
  where
    source = forM_ fls $ \fl -> case fl of
        Left fl'  -> sourceFileBS (fl'^.location)
        Right fl' -> sourceFileBS (fl'^.location) .| ungzip
{-# INLINE concatBed #-}
