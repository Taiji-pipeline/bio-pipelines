{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
module Bio.Pipeline.NGS.STAR
    ( STAROpts
    , defaultSTAROpts
    , starCmd
    , starCores
    , starSort
    , starTranscriptome
    , starTmpDir
    , starMkIndex
    , starAlign
    ) where

import           Bio.Data.Experiment
import           Control.Lens                ((.~), (^.), makeLenses)
import           Control.Monad.State.Lazy
import           Data.Either                 (isRight)
import           Data.Maybe                  (isJust)
import           Data.Singletons.Prelude.List (Delete)
import           Data.Singletons             (SingI)
import qualified Data.Text                   as T
import           Shelly                      (fromText, mkdir_p, cp, run_,
                                              shelly, test_f)
import           System.IO                   (hPutStrLn, stderr)
import           System.IO.Temp              (withTempDirectory)

data STAROpts = STAROpts
    { _starCmd           :: FilePath    -- ^ The path to STAR executable
    , _starCores         :: Int         -- ^ Number of parallel threads
    , _starTmpDir        :: FilePath    -- ^ Directory for storing temporary files
    , _starSort          :: Bool        -- ^ Whether to sort output Bam file
    , _starTranscriptome :: Maybe FilePath -- ^ Output in transcript coordinates
    }

makeLenses ''STAROpts

defaultSTAROpts :: STAROpts
defaultSTAROpts = STAROpts
    { _starCmd = "STAR"
    , _starCores = 1
    , _starTmpDir = "./"
    , _starSort = False
    , _starTranscriptome = Nothing
    }

-- | Create index files for STAR
starMkIndex :: FilePath   -- ^ STAR command path
            -> FilePath   -- ^ Directory used to store genome indices
            -> [FilePath] -- ^ Fastq files
            -> FilePath   -- ^ Annotation file
            -> Int        -- ^ The length of the genomic sequence
                          -- around the annotated junction to be used in
                          -- constructing the splice junctions database. Set it
                          -- to "ReadLength-1" or 100 for general purpose.
            -> IO FilePath
starMkIndex star idx fstqs anno r = do
    indexExist >>= \case
        True -> hPutStrLn stderr "STAR index directory exists. Skipped."
        False -> shelly $ do
            mkdir_p idx
            liftIO $ hPutStrLn stderr "Generating STAR indices"
            run_ star $
                [ "--runThreadN", "4", "--runMode", "genomeGenerate", "--genomeDir"
                , T.pack idx, "--genomeFastaFiles" ] ++ map T.pack fstqs ++
                ["--sjdbGTFfile", T.pack anno, "--sjdbOverhang", T.pack $ show r]
    return idx
  where
    indexExist = shelly $ fmap and $ forM exts $ \ext -> test_f $ idx <> "/" <> ext
      where
        exts = ["chrLength.txt", "chrName.txt", "chrStart.txt"
            , "Genome", "genomeParameters.txt", "SA", "SAindex"]

-- | Align RNA-seq raw reads with STAR
starAlign :: ( SingI tags
             , tags1 ~ Delete 'Gzip tags
             , tags2 ~ Insert' 'PairedEnd tags1 )
          => FilePath                    -- ^ Genome alignment result
          -> FilePath                    -- ^ STAR genome index
          -> Either (File tags 'Fastq)
                    (File tags 'Fastq, File tags 'Fastq)
          -> STAROpts                    -- ^ Options
          -> IO ( Either (File tags1 'Bam, Maybe (File tags1 'Bam))
                         (File tags2 'Bam, Maybe (File tags2 'Bam)) )
starAlign outputGenome index dat opt = withTempDirectory
    (opt^.starTmpDir) "STAR_align_tmp_dir." $ \tmp_dir -> shelly $ do
        run_ (opt^.starCmd) $ ["--readFilesIn"] ++ map T.pack inputs ++
            ["--genomeDir", T.pack index
            , "--outFileNamePrefix", T.pack $ tmp_dir ++ "/"
            , "--runThreadN",  T.pack $ show $ opt^.starCores
            , "--genomeLoad", "NoSharedMemory"
            , "--outFilterType", "BySJout"     -- reduces the number of ”spurious” junctions
            , "--outFilterMultimapNmax", "20"  -- max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
            , "--alignSJoverhangMin", "8"      -- minimum overhang for unannotated junctions
            , "--alignSJDBoverhangMin", "1"    -- minimum overhang for annotated junctions
            , "--outFilterMismatchNmax", "999" -- maximum number of mismatches per pair, large number switches off this filter
            , "--outFilterMismatchNoverReadLmax", "0.04" -- max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.06*200=8 for the paired read
            , "--alignIntronMin", "20"        -- minimum intron length
            , "--alignIntronMax", "1000000"    -- maximum intron length
            , "--alignMatesGapMax", "1000000"  -- maximum genomic distance between mates
            , "--outSAMunmapped", "Within", "--outSAMattributes"
            , "NH", "HI", "AS", "NM", "MD"
            , "--outSAMheaderCommentFile", "COfile.txt"
            , "--outSAMheaderHD", "@HD", "VN:1.4", "SO:coordinate"
            , "--sjdbScore", "1" ] ++
            ( if zipped then ["--readFilesCommand", "zcat"] else [] ) ++
            ( if opt^.starSort
                then [ "--outSAMtype", "BAM", "SortedByCoordinate"
                     , "--limitBAMsortRAM", "60000000000" ]
                else ["--outSAMtype", "BAM", "Unsorted"] ) ++
            ( if isRight dat then [] else ["--outSAMstrandField", "intronMotif"] ) ++
            ( if isJust (opt^.starTranscriptome)
                then ["--quantMode", "TranscriptomeSAM"]
                else [] )

        let starOutput | opt^.starSort = "/Aligned.sortedByCoord.out.bam"
                       | otherwise = "/Aligned.out.bam"
        cp (tmp_dir ++ starOutput) outputGenome

        -- Sorting annotation bam
        case opt^.starTranscriptome of
            Just outputAnno -> do
                if opt^.starSort
                    then run_ "samtools" [ "sort", "-@", T.pack $ show $ opt^.starCores
                            , "-T", T.pack $ tmp_dir ++ "/sort_bam_tmp"
                            , "-o", T.pack outputAnno
                            , T.pack $ tmp_dir ++ "/Aligned.toTranscriptome.out.bam" ]
                    else cp (tmp_dir ++ "/Aligned.toTranscriptome.out.bam") outputAnno
                return $ case dat of
                    Left _ -> Left ( location .~ outputGenome $ emptyFile
                                   , Just $ location .~ outputAnno $ emptyFile )
                    Right _ -> Right ( location .~ outputGenome $ emptyFile
                                     , Just $ location .~ outputAnno $ emptyFile )
            Nothing -> do
                return $ case dat of
                    Left _ -> Left (location .~ outputGenome $ emptyFile, Nothing)
                    Right _ -> Right (location .~ outputGenome $ emptyFile, Nothing)
  where
    (inputs, zipped) = case dat of
        Left fastq     -> ([fastq^.location], fastq `hasTag` Gzip)
        Right (f1, f2) -> ([f1^.location, f2^.location], f1 `hasTag` Gzip)