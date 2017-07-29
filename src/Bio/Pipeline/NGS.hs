{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE ExtendedDefaultRules #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE OverloadedLists      #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE TemplateHaskell      #-}

module Bio.Pipeline.NGS
    ( BWAOpts
    , BWAOptSetter
    , bwaCores
    , bwaSeedLen
    , bwaMaxMis
    , bwaReadTrim
    , defaultBWAOpts
    , bwaMkIndex
    , bwaAlign
    , filterBam
    , removeDuplicates
    , bam2Bed
    , bam2BedPE
    , mergeReplicatesBed

    , STAROpts
    , STAROptSetter
    , starCmd
    , starCores
    , starSort
    , starTmpDir
    , starMkIndex
    , starAlign

    , RSEMOpts
    , RSEMOptSetter
    , rsemPath
    , rsemCores
    , rsemSeed
    , rsemMkIndex
    , rsemQuant
    ) where

import           Bio.Data.Bam             (bamToBed, readBam, runBam,
                                           sortedBamToBedPE)
import           Bio.Data.Bed             (BED, BED3 (..), BEDLike (..), toLine)
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib        (gzip, ungzip)
import           Data.Maybe               (fromJust)
import           Data.Tagged              (Tagged (..), untag)
import qualified Data.Text                as T
import           Shelly                   (escaping, fromText, mv, run_, shelly,
                                           silently)
import           System.IO.Temp           (withTempDirectory)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR

-- | Remove low quality and redundant tags, fill in mate information.
filterBam :: Bool      -- ^ pairedEnd
          -> FilePath  -- ^ output
          -> File 'Bam -> IO (File 'Bam)
filterBam isPair output fl = withTempDirectory "./" "tmp_filt_dir." $ \tmp -> do
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

-- | Remove duplicates
removeDuplicates :: Bool -> FilePath -> FilePath -> File 'Bam
                 -> IO (File 'Bam, File 'Other)
removeDuplicates isPair picardPath output input =
    withTempDirectory "./" "tmp_picard_dir." $ \tmp -> shelly $ do
        let qcFile = output ++ ".picard.qc"
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
        if isPair
            then run_ "samtools" [ "sort", T.pack filtTmp, "-n", "-T"
                , T.pack $ tmp ++ "/tmp_sort", "-o", T.pack output ]
            else mv (fromText $ T.pack filtTmp) $ fromText $ T.pack output

        let finalBam = tags .~ ["processed bam file"] $ location .~ output $ emptyFile
            dupQC =tags .~ ["picard qc file"] $ location .~ qcFile $ emptyFile
        return (finalBam, dupQC)

bam2Bed :: String    -- ^ Prefix
        -> (BED -> Bool)  -- ^ Filtering function
        -> File 'Bam -> IO (Tagged GZipped (File 'Bed))
bam2Bed output fn fl = do
    runBam $ readBam (fl^.location) =$= bamToBed =$= filterC fn =$=
        mapC toLine =$= unlinesAsciiC =$= gzip $$ sinkFileBS output
    return $ Tagged $ location .~ output $ emptyFile
{-# INLINE bam2Bed #-}

-- | Convert name sorted BAM to BEDPE suitable for MACS2.
bam2BedPE :: String
          -> ((BED, BED) -> Bool)
          -> Tagged Sorted (File 'Bam)
          -> IO (Tagged GZipped (File 'Bed))
bam2BedPE output fn fl = do
    runBam $ readBam (untag fl ^. location) =$= sortedBamToBedPE =$=
        filterC fn =$= concatMapC f =$= mapC toLine =$= unlinesAsciiC =$=
        gzip $$ sinkFileBS output
    return $ Tagged $ location .~ output $ emptyFile
  where
    f (b1, b2)
        | chrom b1 /= chrom b2 || bedStrand b1 == bedStrand b2 = Nothing
        | otherwise =
            let left = if fromJust (bedStrand b1)
                    then chromStart b1 else chromStart b2
                right = if not (fromJust $ bedStrand b2)
                    then chromEnd b2 else chromEnd b1
            in if left < right
                  then Just $ BED3 (chrom b1) left right
                  else error "Left coordinate is larger than right coordinate."
{-# INLINE bam2BedPE #-}

-- | Merge multiple BED files.
mergeReplicatesBed :: FilePath -> [MaybeTagged GZipped (File 'Bed)]
                   -> IO (Tagged GZipped (File 'Bed))
mergeReplicatesBed output fls = do
    let source = forM_ fls $ \fl -> case fl of
            Right (Tagged x) -> sourceFileBS (x^.location) =$= ungzip
            Left x           -> sourceFileBS (x^.location)
    runResourceT $ source =$= gzip $$ sinkFile output
    return $ Tagged $ location .~ output $ emptyFile
{-# INLINE mergeReplicatesBed #-}
