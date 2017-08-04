{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Bio.Pipeline.NGS.STAR
    ( STAROpts
    , STAROptSetter
    , defaultSTAROpts
    , starCmd
    , starCores
    , starSort
    , starTmpDir
    , starMkIndex
    , starAlign
    ) where

import           Bio.Data.Experiment
import           Control.Lens             ((.~), (^.))
import           Control.Lens             (makeLenses)
import           Control.Monad.State.Lazy
import qualified Data.Text                as T
import           Shelly                   (fromText, mkdir_p, mv, run_, shelly,
                                           test_d)
import           System.IO                (hPutStrLn, stderr)
import           System.IO.Temp           (withTempDirectory)

data STAROpts = STAROpts
    { _starCmd    :: FilePath
    , _starCores  :: Int
    , _starTmpDir :: FilePath
    , _starSort   :: Bool
    }

makeLenses ''STAROpts

type STAROptSetter = State STAROpts ()

defaultSTAROpts :: STAROpts
defaultSTAROpts = STAROpts
    { _starCmd = "STAR"
    , _starCores = 1
    , _starTmpDir = "./"
    , _starSort = False
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
starMkIndex star dir fstqs anno r = do
    dirExist <- shelly $ test_d $ fromText $ T.pack dir
    if dirExist
        then hPutStrLn stderr "STAR index directory exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack dir
            liftIO $ hPutStrLn stderr "Generating STAR indices"
            run_ (fromText $ T.pack star) $
                [ "--runThreadN", "1", "--runMode", "genomeGenerate", "--genomeDir"
                , T.pack dir, "--genomeFastaFiles" ] ++ map T.pack fstqs ++
                ["--sjdbGTFfile", T.pack anno, "--sjdbOverhang", T.pack $ show r]
    return dir

-- | Align RNA-seq raw reads with STAR
starAlign :: ( MayHave 'Gzip tags
             , MayHave Pairend tags
             , tags' ~ Remove 'Gzip tags )
          => FilePath                    -- ^ Genome alignment result
          -> FilePath                    -- ^ Annotation alignment result
          -> FilePath                    -- ^ STAR genome index
          -> STAROptSetter               -- ^ Options
          -> Either (File tags 'Fastq)
                    (File tags 'Fastq, File tags 'Fastq)
          -> IO (File tags' 'Bam, File tags' 'Bam)
starAlign outputGenome outputAnno index setter dat = withTempDirectory
    (opt^.starTmpDir) "STAR_align_tmp_dir." $ \tmp_dir -> shelly $ do
        run_ star $ ["--genomeDir", T.pack index, "--readFilesIn"] ++
            map (T.pack . (^.location)) inputs ++
            [ "--outFileNamePrefix", T.pack $ tmp_dir ++ "/"
            , "--runThreadN",  T.pack $ show $ opt^.starCores ] ++
            ( if zipped then ["--readFilesCommand", "zcat"] else [] ) ++
            [ "--genomeLoad", "NoSharedMemory"
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
            , "--outSAMheaderHD", "@HD", "VN:1.4", "SO:coordinate" ] ++
            ( if opt^.starSort
                then [ "--outSAMtype", "BAM", "SortedByCoordinate"
                     , "--limitBAMsortRAM", "60000000000" ]
                else ["--outSAMtype", "BAM", "Unsorted"] ) ++
            ( if isPair then [] else ["--outSAMstrandField", "intronMotif"] ) ++
            ["--quantMode", "TranscriptomeSAM", "--sjdbScore", "1"]

        let starOutput | opt^.starSort = "/Aligned.sortedByCoord.out.bam"
                       | otherwise = "/Aligned.out.bam"
        mv (fromText $ T.pack $ tmp_dir ++ starOutput) $ fromText $
            T.pack $ outputGenome

        -- Sorting annotation bam
        if opt^.starSort
            then run_ "samtools" [ "sort", "-@", T.pack $ show $ opt^.starCores
                    , "-T", T.pack $ tmp_dir ++ "/sort_bam_tmp"
                    , "-o", T.pack outputAnno
                    , T.pack $ tmp_dir ++ "/Aligned.toTranscriptome.out.bam" ]
            else mv (fromText $ T.pack $ tmp_dir ++
                    "/Aligned.toTranscriptome.out.bam") $ fromText $
                    T.pack outputAnno

        let genomeAlignFile = location .~ outputGenome $
                tags .~ ["RNA genome align bam"] $ emptyFile
            annoFile = location .~ outputAnno $
                tags .~ ["RNA anno align bam"] $ emptyFile
        return (genomeAlignFile, annoFile)
  where
    star = fromText $ T.pack $ opt^.starCmd
    (inputs, zipped, isPair) = case dat of
        Left fastq -> ( [fastq], isGzipped fastq, False)
        Right (f1, f2) -> ([f1,f2], isGzipped f1
            , if isPairend f1 then True else error "Must be pairended")
    opt = execState setter defaultSTAROpts
