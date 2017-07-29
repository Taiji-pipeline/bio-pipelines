{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedLists        #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}
{-# LANGUAGE DataKinds #-}

module Bio.Pipeline.CallPeaks
    ( CallPeakOpts(..)
    , CallPeakOptSetter
    , Cutoff(..)
    , tmpDir
    , cutoff
    , gSize
    , pair
    , defaultCallPeakOpts
    , callPeaks
    , idr
    , idrMultiple
    ) where

import qualified Bio.Data.Bed as Bed
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import qualified Data.ByteString.Char8     as B
import           Data.Conduit.Zlib         (ungzip)
import           Data.Ord
import qualified Data.Text                 as T
import           Shelly                    (fromText, mv, run_, shelly)
import           System.IO.Temp            (withTempDirectory)

import           Data.List

type CallPeakOptSetter = State CallPeakOpts ()

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir :: FilePath
    , callPeakOptsCutoff :: Cutoff
    , callPeakOptsGSize  :: String
    , callPeakOptsPair   :: Bool
    --, callPeakOptsBroad :: Bool
    --, callPeakOptsBroadCutoff :: Double
    }

data Cutoff = PValue Double
            | QValue Double

makeFields ''CallPeakOpts

defaultCallPeakOpts :: CallPeakOpts
defaultCallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir = "./"
    , callPeakOptsCutoff = QValue 0.01
    , callPeakOptsGSize = "mm"
    , callPeakOptsPair  = False
    --, callPeakOptsBroad = False
    --, callPeakOptsBroadCutoff = 0.05
    }

-- | Call peaks using MACS2.
callPeaks :: FilePath           -- ^ Ouptut file
          -> File 'Bed    -- ^ Sample
          -> Maybe (File 'Bed)      -- ^ Input/control sample
          -> CallPeakOptSetter  -- ^ Options
          -> IO (File 'NarrowPeak)
callPeaks output target input setter = do
    macs2 output (target^.location) (fmap (^.location) input)
        fileFormat opt
    f <- frip (target^.location) output
    return $ location .~ output $
        tags .~ ["macs2"] $ info .~ [("FRiP", T.pack $ show f)] $ emptyFile
  where
    opt = execState setter defaultCallPeakOpts
    fileFormat | opt^.pair = "BEDPE"
               | otherwise = "AUTO"
{-# INLINE callPeaks #-}

macs2 :: FilePath        -- ^ Output
      -> FilePath        -- ^ Target
      -> Maybe FilePath  -- ^ Input
      -> String          -- ^ File format
      -> CallPeakOpts
      -> IO ()
macs2 output target input fileformat opt = withTempDirectory (opt^.tmpDir)
    "tmp_macs2_dir." $ \tmp -> shelly $ do
        run_ "macs2" $
            [ "callpeak", "-f", T.pack fileformat, "-g", T.pack $ opt^.gSize
            , "--outdir", T.pack tmp, "--tempdir", T.pack tmp, "--keep-dup"
            , "all", "-t", T.pack target
            ] ++ control ++ cut
        mv (fromText $ T.pack $ tmp ++ "/NA_peaks.narrowPeak") $ fromText $
            T.pack output
  where
    control = case input of
        Nothing -> []
        Just x -> ["-c", T.pack x]
    cut = case opt^.cutoff of
        QValue x -> ["--qvalue", T.pack $ show x]
        PValue x -> ["--pvalue", T.pack $ show x]
{-# INLINE macs2 #-}

-- | Fraction of reads in peaks
frip :: FilePath   -- ^ reads, in BedGzip format
     -> FilePath   -- ^ peaks, in bed format
     -> IO Double
frip rs peaks = do
    p <- Bed.readBed' peaks :: IO [Bed.BED3]
    (n, m) <- flip execStateT (0::Int, 0::Int) $ runResourceT $
        sourceFileBS rs =$= ungzip =$= linesUnboundedAsciiC =$=
        mapC (Bed.fromLine :: B.ByteString -> Bed.BED3) =$= total =$=
        Bed.intersectBed p $$ count
    return $ fromIntegral m / fromIntegral n
  where
    total = awaitForever $ \i -> do
        (c, x) <- get
        put (c+1,x)
        yield i
    count = awaitForever $ \_ -> do
        (x, c) <- get
        put (x, c+1)

idrMultiple :: [File 'NarrowPeak]   -- ^ Peaks
            -> File 'NarrowPeak  -- ^ Merged peaks
            -> Double
            -> FilePath
            -> IO (File 'NarrowPeak)
idrMultiple [x] _ _ _ = return x
idrMultiple peakFiles merged th output =
    withTempDirectory "./" "tmp_idr_dir." $ \tmp -> do
        peaks <- forM (zip [1::Int ..] peakPair) $ \(i, (p1, p2)) -> do
            result <- idr p1 p2 merged th $ tmp ++ "/" ++ show i
            n <- numLine $ result^.location
            return (result^.location, n)
        let final = fst $ maximumBy (comparing snd) peaks
        shelly $ mv (fromText $ T.pack final) $ fromText $ T.pack output
        return $ location .~ output $ tags .~ ["IDR"] $ emptyFile
  where
    numLine x = do
        c <- B.readFile x
        return $ length $ B.lines c
    peakPair = comb peakFiles
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _ = []

-- | Perform Irreproducible Discovery Rate (IDR) analysis
idr :: File 'NarrowPeak  -- ^ Peak 1
    -> File 'NarrowPeak  -- ^ Peak 2
    -> File 'NarrowPeak  -- ^ Peaks called from merged replicates (relax threshold)
    -> Double    -- ^ IDR threshold
    -> FilePath  -- ^ Output
    -> IO (File 'NarrowPeak)
idr peak1 peak2 peakMerged th output = do
    shelly $ run_ "idr" [ "--samples", p1, p2, "--peak-list", pm
        , "--input-file-type", "narrowPeak", "--rank", "signal.value"
        , "--idr-threshold", T.pack $ show th, "-o", T.pack output ]
    return $ location .~ output $
        tags .~ ["IDR"] $ emptyFile
  where
    p1 = T.pack $ peak1^.location
    p2 = T.pack $ peak2^.location
    pm = T.pack $ peakMerged^.location
