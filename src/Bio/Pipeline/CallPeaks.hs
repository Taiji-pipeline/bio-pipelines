{-# LANGUAGE DataKinds              #-}
{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedLists        #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}
{-# LANGUAGE LambdaCase #-}

module Bio.Pipeline.CallPeaks
    ( CallPeakOpts(..)
    , CallPeakMode(..)
    , Cutoff(..)
    , tmpDir
    , cutoff
    , gSize
    , mode
    , callSummits
    , genSignal
    , noLambda
    , rawParam
    , callPeaks
    , frip
    , idr
    , idrMultiple
    ) where

import qualified Bio.Data.Bed          as Bed
import           Bio.Data.Experiment
import           Conduit
import Data.Conduit.Internal (zipSinks)
import           Control.Lens
import           Control.Monad
import qualified Data.ByteString.Char8 as B
import           Data.Conduit.Zlib     (gzip)
import           Data.Default          (Default (..))
import           Data.List
import           Data.Ord
import           Data.Singletons       (SingI)
import qualified Data.Text             as T
import           Shelly                (fromText, cp, run_, shelly)
import           System.IO.Temp        (withTempDirectory)

data CallPeakMode = Model
                  | NoModel Int Int   -- ^ implies "--nomodel --shift n --extsize m"

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir      :: FilePath
    , callPeakOptsCutoff      :: Cutoff
    , callPeakOptsGSize       :: Maybe String
    , callPeakOptsMode        :: CallPeakMode
    , callPeakOptsCallSummits :: Bool
    , callPeakOptsGenSignal   :: Bool
    , callPeakOptsNoLambda    :: Bool
    , callPeakOptsRawParam :: [T.Text]
    --, callPeakOptsBroad :: Bool
    --, callPeakOptsBroadCutoff :: Double
    }

data Cutoff = PValue Double
            | QValue Double

makeFields ''CallPeakOpts

instance Default CallPeakOpts where
    def = CallPeakOpts
        { callPeakOptsTmpDir = "./"
        , callPeakOptsCutoff = QValue 0.01
        , callPeakOptsGSize = Nothing
        , callPeakOptsMode  = Model
        , callPeakOptsCallSummits = False
        , callPeakOptsGenSignal = False
        , callPeakOptsNoLambda = False
        , callPeakOptsRawParam = []
        --, callPeakOptsBroad = False
        --, callPeakOptsBroadCutoff = 0.05
        }

-- | Call peaks using MACS2.
callPeaks :: SingI tags
          => FilePath                  -- ^ Ouptut file
          -> File tags 'Bed            -- ^ Sample
          -> Maybe (File tags 'Bed)    -- ^ Input/control sample
          -> CallPeakOpts              -- ^ Options
          -> IO (File '[Gzip] 'NarrowPeak)
callPeaks output target input opt = runResourceT
    (runConduit $ streamer .| nullC) >>= \case
        True -> error "Call peaks failed because the input BED file is empty"
        False -> do
            macs2 output (target^.location) (fmap (^.location) input)
                fileFormat opt
            return $ location .~ output $ emptyFile
  where
    fileFormat | target `hasTag` PairedEnd = "BEDPE"
               | otherwise = "BED"
    streamer :: ConduitT () Bed.BED3 (ResourceT IO) ()
    streamer = if target `hasTag` Gzip
        then Bed.streamBedGzip (target^.location)
        else Bed.streamBed (target^.location)

macs2 :: FilePath        -- ^ Output
      -> FilePath        -- ^ Target
      -> Maybe FilePath  -- ^ Input
      -> String          -- ^ File format
      -> CallPeakOpts
      -> IO ()
macs2 output target input fileformat opt = withTempDirectory (opt^.tmpDir)
    "tmp_macs2_dir." $ \tmp -> shelly $ do
        run_ "macs2" $ [ "callpeak"
            , "-f", T.pack fileformat
            , "--outdir", T.pack tmp
            , "--tempdir", T.pack tmp
            , "--keep-dup", "all"
            , "-t", T.pack target ]
            ++ ( case opt^.gSize of
                    Nothing -> []
                    Just x ->  ["-g", T.pack x] )
            ++ ( case input of
                    Nothing -> []
                    Just x  -> ["-c", T.pack x] )
            ++ (if opt^.callSummits then ["--call-summits"] else [])
            ++ ( case opt^.cutoff of
                    QValue x -> ["--qvalue", T.pack $ show x]
                    PValue x -> ["--pvalue", T.pack $ show x] )
            ++ ( case opt^.mode of
                    Model -> []
                    NoModel shift ext ->
                        [ "--nomodel", "--shift", T.pack $ show shift
                        , "--extsize", T.pack $ show ext ] )
            ++ (if opt^.genSignal then ["-B", "--SPMR"] else [])
            ++ (if opt^.noLambda then ["--nolambda"] else [])
            ++ (opt^.rawParam)
        liftIO $ runResourceT $ runConduit $
            sourceFileBS (tmp <> "/NA_peaks.narrowPeak") .| gzip .| sinkFile output
        when (opt^.genSignal) $ do
            let treatBdg = T.pack $ tmp ++ "/NA_treat_pileup.bdg"
                ctrlBdg = T.pack $ tmp ++ "/NA_control_lambda.bdg"
            shelly $ run_ "macs2" ["bdgcmp", "-t", treatBdg, "-c", ctrlBdg,
                "-m", "FE", "-o", T.pack $ output <> ".bdg"]
{-# INLINE macs2 #-}

-- | Fraction of reads in peaks
frip :: (SingI tags1, SingI tags2)
     => File tags1 'Bed        -- ^ reads
     -> File tags2 'NarrowPeak -- ^ peaks
     -> IO Double
frip rs peak = do
    p <- runResourceT $ runConduit $ streamer peak .| sinkList
    (m, n) <- runResourceT $ runConduit $ streamer rs .| zipSinks (Bed.intersectBed p .| lengthC) lengthC
    return $ fromIntegral (m :: Int) / fromIntegral (n :: Int)
  where
    streamer :: SingI tags => File tags file -> ConduitT () Bed.BED3 (ResourceT IO) ()
    streamer x = if x `hasTag` Gzip then Bed.streamBedGzip (x^.location) else Bed.streamBed (x^.location)

idrMultiple :: [File tags 'NarrowPeak]   -- ^ Peaks
            -> File tags 'NarrowPeak  -- ^ Merged peaks
            -> Double
            -> FilePath
            -> IO (File tags 'NarrowPeak)
idrMultiple [x] _ _ _ = return x
idrMultiple peakFiles merged th output =
    withTempDirectory "./" "tmp_idr_dir." $ \tmp -> do
        peaks <- forM (zip [1::Int ..] peakPair) $ \(i, (p1, p2)) -> do
            result <- idr p1 p2 merged th $ tmp ++ "/" ++ show i
            n <- numLine $ result^.location
            return (result^.location, n)
        let final = fst $ maximumBy (comparing snd) peaks
        shelly $ cp (fromText $ T.pack final) $ fromText $ T.pack output
        return $ location .~ output $ emptyFile
  where
    numLine x = do
        c <- B.readFile x
        return $ length $ B.lines c
    peakPair = comb peakFiles
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _      = []

-- | Perform Irreproducible Discovery Rate (IDR) analysis
idr :: File tags 'NarrowPeak  -- ^ Peak 1
    -> File tags 'NarrowPeak  -- ^ Peak 2
    -> File tags 'NarrowPeak  -- ^ Peaks called from merged replicates (relax threshold)
    -> Double    -- ^ IDR threshold
    -> FilePath  -- ^ Output
    -> IO (File tags 'NarrowPeak)
idr peak1 peak2 peakMerged th output = do
    shelly $ run_ "idr" [ "--samples", p1, p2, "--peak-list", pm
        , "--input-file-type", "narrowPeak", "--rank", "signal.value"
        , "--idr-threshold", T.pack $ show th, "-o", T.pack output ]
    return $ location .~ output $ emptyFile
  where
    p1 = T.pack $ peak1^.location
    p2 = T.pack $ peak2^.location
    pm = T.pack $ peakMerged^.location
