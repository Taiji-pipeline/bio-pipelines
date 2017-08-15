{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Bio.Pipeline.NGS.RSEM
    ( RSEMOpts
    , RSEMOptSetter
    , rsemPath
    , rsemCores
    , rsemSeed
    , rsemMkIndex
    , rsemQuant_
    ) where

import           Bio.Data.Experiment
import           Control.Lens             ((.~), (^.))
import           Control.Lens             (makeLenses)
import           Control.Monad.State.Lazy
import           Data.Int                 (Int32)
import           Data.Singletons          (SingI)
import qualified Data.Text                as T
import           Shelly                   (fromText, mkdir_p, run_, shelly,
                                           test_d)
import           System.FilePath          (takeDirectory)
import           System.IO                (hPutStrLn, stderr)

rsemMkIndex :: FilePath   -- ^ Prefix
            -> FilePath   -- ^ annotation file in GFF3 format
            -> [FilePath] -- ^ fastq files
            -> IO FilePath
rsemMkIndex prefix anno fstqs = do
    dirExist <- shelly $ test_d $ fromText $ T.pack dir
    if dirExist
        then hPutStrLn stderr "RSEM index directory exists. Skipped."
        else shelly $ do
            mkdir_p $ fromText $ T.pack dir
            liftIO $ hPutStrLn stderr "Generating RSEM indices"
            run_ "rsem-prepare-reference" [ "--gtf", T.pack anno
                , T.intercalate "," $ map T.pack fstqs, T.pack prefix ]
    return prefix
  where
    dir = takeDirectory prefix


data RSEMOpts = RSEMOpts
    { _rsemPath  :: FilePath
    , _rsemCores :: Int
    , _rsemSeed  :: Int32
    }

makeLenses ''RSEMOpts

type RSEMOptSetter = State RSEMOpts ()

defaultRSEMOpts :: RSEMOpts
defaultRSEMOpts = RSEMOpts
    { _rsemPath = ""
    , _rsemCores = 1
    , _rsemSeed = 12345
    }

-- | Gene and transcript quantification using rsem
rsemQuant_ :: SingI tags
           => FilePath         -- ^ output prefix
           -> FilePath         -- ^ Directory containing the index
           -> RSEMOptSetter
           -> File tags 'Bam
           -> IO (File tags 'Tsv, File tags 'Tsv)
rsemQuant_ outputPrefix indexPrefix setter input = shelly $ do
    run_ rsem $ [ "--bam", "--estimate-rspd", "--calc-ci"
        , "--seed", T.pack $ show $ opt^.rsemSeed
        , "-p", T.pack $ show $ opt^.rsemCores
        , "--no-bam-output", "--ci-memory", "30000" ] ++
        ( if input `hasTag` Pairend
            then ["--paired-end", "--forward-prob", "0"]
            else [] ) ++
        [T.pack $ input^.location, T.pack indexPrefix, T.pack outputPrefix]

    let geneQuant = location .~ outputPrefix ++ ".genes.results" $ emptyFile
        transcirptQuant = location .~ outputPrefix ++ ".isoforms.results" $ emptyFile
    return (geneQuant, transcirptQuant)
  where
    rsem = fromText $ T.pack $ opt^.rsemPath ++ "rsem-calculate-expression"
    opt = execState setter defaultRSEMOpts
