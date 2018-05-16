{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Bio.Pipeline.NGS.BWA
    ( BWAOpts
    , defaultBWAOpts
    , bwaCores
    , bwaSeedLen
    , bwaTmpDir
    , bwaMkIndex
    , bwaAlign
    ) where

import           Bio.Data.Experiment
import           Control.Lens                ((.~), (^.))
import           Control.Lens                (makeLenses)
import           Control.Monad.State.Lazy
import           Data.Promotion.Prelude.List (Delete)
import qualified Data.Text                   as T
import           Shelly                      (cp, escaping, fromText, mkdir_p,
                                              run_, shelly, test_f, touchfile)
import           System.FilePath             (takeDirectory)
import           System.IO                   (hPutStrLn, stderr)

data BWAOpts = BWAOpts
    { _bwaCores   :: Int       -- ^ number of cpu cores
    , _bwaSeedLen :: Int       -- ^ seed length, equivalent to -k
    , _bwaTmpDir  :: FilePath  -- ^ temp dir
    } deriving (Show)

makeLenses ''BWAOpts

defaultBWAOpts :: BWAOpts
defaultBWAOpts = BWAOpts
    { _bwaCores = 1
    , _bwaSeedLen = 32
    , _bwaTmpDir = "./"
    }

-- | Generate BWA genome index
bwaMkIndex :: FilePath
           -> FilePath   -- ^ Index prefix, e.g., /path/genome.fa
           -> IO FilePath
bwaMkIndex input prefix = do
    fileExist <- shelly $ test_f $ fromText $ T.pack $ dir ++ stamp
    if fileExist
        then hPutStrLn stderr "BWA index exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack dir
            cp (fromText $ T.pack input) $ fromText $ T.pack prefix
            liftIO $ hPutStrLn stderr "Generating BWA index"
            run_ "bwa" ["index", "-p", T.pack prefix, "-a", "bwtsw", T.pack input]
            touchfile $ fromText $ T.pack $ dir ++ stamp
    return prefix
  where
    dir = takeDirectory prefix
    stamp = "/.bio_pipelines_bwa_index"

bwaAlign :: FilePath  -- ^ Path for the output bam file
         -> FilePath  -- ^ Genome index
         -> Either (File tags 'Fastq)
                   (File tags 'Fastq, File tags 'Fastq)
         -> BWAOpts
         -> IO ( Either (File (Delete 'Gzip tags) 'Bam)
                        (File (Insert' 'PairedEnd (Delete 'Gzip tags)) 'Bam) )
bwaAlign output index fastq opt = do
    shelly $ escaping False $ run_ "bwa" $
        [ "mem", "-M"  -- "picard compatibility"
        , "-k", T.pack $ show $ opt^.bwaSeedLen
        , "-t", T.pack $ show $ opt^.bwaCores
        , T.pack index ] ++ inputs ++
        [ "|", "samtools", "view", "-Su", "-"
        , ">", T.pack $ output ]
    return $ case fastq of
        Left _ -> Left $ location .~ output $ emptyFile
        Right _ -> Right $ location .~ output $ emptyFile
  where
    inputs = case fastq of
        Left f        -> [T.pack $ f^.location]
        Right (f1,f2) -> [T.pack $ f1^.location, T.pack $ f2^.location]
{-# INLINE bwaAlign #-}
