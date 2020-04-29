{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Bio.Pipeline.NGS.RSEM
    ( RSEMOpts
    , defaultRSEMOpts
    , rsemPath
    , rsemCores
    , rsemSeed
    , rsemMkIndex
    , rsemQuant
    ) where

import           Bio.Data.Experiment
import           Control.Lens             ((.~), (^.))
import           Control.Lens             (makeLenses)
import           Control.Monad.State.Lazy
import           Data.Int                 (Int32)
import           Data.Singletons          (SingI)
import qualified Data.Text                as T
import           Shelly                   (fromText, mkdir_p, run_, shelly,
                                           test_f)
import           System.FilePath          (takeDirectory)
import           System.IO                (hPutStrLn, stderr)

rsemMkIndex :: FilePath   -- ^ Prefix
            -> FilePath   -- ^ annotation file in GFF3 format
            -> [FilePath] -- ^ fastq files
            -> IO FilePath
rsemMkIndex prefix anno fstqs = do
    indexExist >>= \case
        True -> hPutStrLn stderr "RSEM index directory exists. Skipped."
        False -> shelly $ do
            mkdir_p $ fromText $ T.pack dir
            liftIO $ hPutStrLn stderr "Generating RSEM indices"
            run_ "rsem-prepare-reference" [ "--gtf", T.pack anno
                , T.intercalate "," $ map T.pack fstqs, T.pack prefix ]
    return prefix
  where
    dir = takeDirectory prefix
    indexExist = shelly $ fmap and $ forM exts $ \ext ->
        test_f $ fromText $ T.pack $ prefix <> ext
      where
        exts = [".idx.fa", ".seq"]

data RSEMOpts = RSEMOpts
    { _rsemPath  :: FilePath
    , _rsemCores :: Int
    , _rsemSeed  :: Int32
    }

makeLenses ''RSEMOpts

defaultRSEMOpts :: RSEMOpts
defaultRSEMOpts = RSEMOpts
    { _rsemPath = ""
    , _rsemCores = 1
    , _rsemSeed = 12345
    }

-- | Gene and transcript quantification using rsem
rsemQuant :: SingI tags
          => FilePath         -- ^ output prefix
          -> FilePath         -- ^ Directory containing the index
          -> File tags 'Bam
          -> RSEMOpts
          -> IO (File '[GeneQuant] 'Tsv, File '[TranscriptQuant] 'Tsv)
rsemQuant outputPrefix indexPrefix input opt = shelly $ do
    run_ rsem $ [ "--bam", "--estimate-rspd", "--calc-ci"
        , "--seed", T.pack $ show $ opt^.rsemSeed
        , "-p", T.pack $ show $ opt^.rsemCores
        , "--no-bam-output", "--ci-memory", "30000" ] ++
        ( if input `hasTag` PairedEnd
            then ["--paired-end", "--forward-prob", "0"]
            else [] ) ++
        [T.pack $ input^.location, T.pack indexPrefix, T.pack outputPrefix]

    let geneQuant = location .~ outputPrefix ++ ".genes.results" $ emptyFile
        transcriptQuant = location .~ outputPrefix ++ ".isoforms.results" $ emptyFile
    return (geneQuant, transcriptQuant)
  where
    rsem = fromText $ T.pack $ opt^.rsemPath ++ "rsem-calculate-expression"
