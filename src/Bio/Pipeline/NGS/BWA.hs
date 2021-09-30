{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE LambdaCase #-}
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
import           Data.Singletons.Prelude.List (Delete)
import qualified Data.Text                   as T
import           Shelly                      (cp, escaping, fromText, mkdir_p,
                                              run_, shelly, test_f,
                                              bashPipeFail, bash_)
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
bwaMkIndex input idx = do
    indexExist >>= \case
        True -> hPutStrLn stderr "BWA index exists. Skipped."
        False -> shelly $ do
            mkdir_p dir
            liftIO $ hPutStrLn stderr "Generating BWA index"
            cp input idx
            run_ "bwa-mem2" ["index", "-p", T.pack idx, T.pack input]
    return idx
  where
    dir = takeDirectory idx
    indexExist = shelly $ fmap and $ forM exts $ \ext -> test_f $ idx <> ext
      where
        exts = ["", ".amb", ".ann", ".pac", ".bwt.2bit.64"]

bwaAlign :: FilePath  -- ^ Path for the output bam file
         -> FilePath  -- ^ Genome index
         -> Either (File tags 'Fastq)
                   (File tags 'Fastq, File tags 'Fastq)
         -> BWAOpts
         -> IO ( Either (File (Delete 'Gzip tags) 'Bam)
                        (File (Insert' 'PairedEnd (Delete 'Gzip tags)) 'Bam) )
bwaAlign output index fastq opt = do
    shelly $ escaping False $ bashPipeFail bash_ "bwa-mem2" $
        [ "mem", "-M"  -- "picard compatibility"
        , "-k", T.pack $ show $ opt^.bwaSeedLen
        , "-t", nCore
        , T.pack index ] ++ inputs ++
        [ "|", "samtools", "view", "--threads", nCore, "-Sb", "-"
        , ">", T.pack $ output ]
    return $ case fastq of
        Left _ -> Left $ location .~ output $ emptyFile
        Right _ -> Right $ location .~ output $ emptyFile
  where
    nCore = T.pack $ show $ opt^.bwaCores
    inputs = case fastq of
        Left f        -> [T.pack $ f^.location]
        Right (f1,f2) -> [T.pack $ f1^.location, T.pack $ f2^.location]
{-# INLINE bwaAlign #-}