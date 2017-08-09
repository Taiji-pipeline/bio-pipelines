{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Bio.Pipeline.NGS.BWA
    ( BWAOpts
    , BWAOptSetter
    , defaultBWAOpts
    , bwaCores
    , bwaSeedLen
    , bwaMaxMis
    , bwaReadTrim
    , bwaTmpDir
    , bwaMkIndex
    , bwaAlign1_
    , bwaAlign2_
    ) where

import           Bio.Data.Experiment
import           Control.Lens                ((.~), (^.))
import           Control.Lens                (makeLenses)
import           Control.Monad.State.Lazy
import           Data.Promotion.Prelude.List (Delete)
import qualified Data.Text                   as T
import           Shelly                      (cp, escaping, fromText, mkdir_p,
                                              run_, shelly, test_f)
import           System.FilePath             (takeDirectory)
import           System.IO                   (hPutStrLn, stderr)
import           System.IO.Temp              (withTempDirectory)

data BWAOpts = BWAOpts
    { _bwaCores    :: Int          -- ^ number of cpu cores
    , _bwaSeedLen  :: Int     -- ^ seed length, equivalent to -l
    , _bwaMaxMis   :: Int    -- ^ max mismatches in seed, equivalent to -k
    , _bwaReadTrim :: Int       -- ^ dynamic read trimming, equivalent to -q
    , _bwaTmpDir   :: FilePath  -- ^ temp dir
    } deriving (Show)

makeLenses ''BWAOpts

defaultBWAOpts :: BWAOpts
defaultBWAOpts = BWAOpts
    { _bwaCores = 1
    , _bwaSeedLen = 32
    , _bwaMaxMis = 2
    , _bwaReadTrim = 5
    , _bwaTmpDir = "./"
    }

type BWAOptSetter = State BWAOpts ()

-- | Generate BWA genome index
bwaMkIndex :: FilePath
           -> FilePath   -- ^ Index prefix, e.g., /path/genome.fa
           -> IO FilePath
bwaMkIndex input prefix = do
    fileExist <- shelly $ test_f (fromText $ T.pack prefix)
    if fileExist
        then hPutStrLn stderr "BWA index exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack $ takeDirectory prefix
            cp (fromText $ T.pack input) $ fromText $ T.pack prefix
            liftIO $ hPutStrLn stderr "Generating BWA index"
            run_ "bwa" ["index", "-p", T.pack prefix, "-a", "bwtsw", T.pack input]
    return prefix

    {-
-- | Tag alignment with BWA aligner.
bwaAlign_ :: (tags1' ~ Delete 'Gzip tags1, tags2' ~ Delete 'Gzip tags2)
          => FilePath  -- ^ Output bam filename
          -> FilePath  -- ^ Genome index
          -> BWAOptSetter
          -> MaybePair tags1 tags2 'Fastq
          -> IO (EitherTag tags1' tags2' 'Bam)
bwaAlign_ output index setter fileset = case fileset of
    Left input             -> Left <$> bwaAlign1_ output index opt input
    Right (input1, input2) -> Right <$> bwaAlign2_ output index opt input1 input2
  where
    opt = execState setter defaultBWAOpts
    -}

bwaAlign1_ :: FilePath  -- ^ Path for the output bam file
           -> FilePath  -- ^ Genome index
           -> BWAOptSetter
           -> File tags 'Fastq
           -> IO (File (Delete 'Gzip tags) 'Bam)
bwaAlign1_ output index setter fastq = do
    stats <- withTempDirectory (opt^.bwaTmpDir) "bwa_align_tmp_dir." $
        \tmpdir -> shelly $ escaping False $ do
            let tmp_sai = T.pack $ tmpdir ++ "/tmp.sai"
            -- Align reads and save the results to tmp_sai.
            run_ "bwa"
                [ "aln", "-q", T.pack $ show $ opt^.bwaReadTrim
                , "-l", T.pack $ show $ opt^.bwaSeedLen
                , "-k", T.pack $ show $ opt^.bwaMaxMis
                , "-t", T.pack $ show $ opt^.bwaCores
                , T.pack index, input, ">", tmp_sai ]
            -- Convert sai to bam.
            run_ "bwa" [ "samse",  T.pack index, tmp_sai, input, "|"
                       , "samtools", "view", "-Su", "-", ">"
                       , T.pack output ]
    return $ location .~ output $ emptyFile
  where
    input = T.pack $ fastq^.location
    opt = execState setter defaultBWAOpts
{-# INLINE bwaAlign1_ #-}

bwaAlign2_ :: FilePath  -- ^ Path for the output bam file
           -> FilePath  -- ^ Genome index
           -> BWAOptSetter
           -> File tags 'Fastq
           -> File tags 'Fastq
           -> IO (File (Delete 'Gzip tags) 'Bam)
bwaAlign2_ output index setter fastqF fastqR = do
    stats <- shelly $ escaping False $
        run_ "bwa" [ "mem", "-M", "-k", T.pack $ show $ opt^.bwaSeedLen
                   , "-t", T.pack $ show $ opt^.bwaCores, T.pack index
                   , input1, input2, "|"
                   , "samtools", "view", "-Su", "-", ">", T.pack $ output ]
    return $ location .~ output $ emptyFile
  where
    opt = execState setter defaultBWAOpts
    input1 = T.pack $ fastqF^.location
    input2 = T.pack $ fastqR^.location
{-# INLINE bwaAlign2_ #-}
