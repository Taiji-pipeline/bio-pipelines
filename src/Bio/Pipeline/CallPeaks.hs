{-# LANGUAGE DataKinds              #-}
{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedLists        #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}

module Bio.Pipeline.CallPeaks
    ( CallPeakOpts(..)
    , CallPeakOptSetter
    , CallPeakMode(..)
    , Cutoff(..)
    , tmpDir
    , cutoff
    , gSize
    , mode
    , callSummits
    , defaultCallPeakOpts
    , callPeaks
    , idr
    , idrMultiple
    ) where

import qualified Bio.Data.Bed             as Bed
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import qualified Data.ByteString.Char8    as B
import           Data.Conduit.Zlib        (ungzip)
import           Data.List
import           Data.Ord
import           Data.Singletons          (SingI)
import qualified Data.Text                as T
import           Shelly                   (fromText, mv, run_, shelly)
import           System.IO.Temp           (withTempDirectory)

data CallPeakMode = Model
                  | NoModel Int Int   -- ^ implies "--nomodel --shift n --extsize m"

type CallPeakOptSetter = State CallPeakOpts ()

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir      :: FilePath
    , callPeakOptsCutoff      :: Cutoff
    , callPeakOptsGSize       :: String
    , callPeakOptsMode        :: CallPeakMode
    , callPeakOptsCallSummits :: Bool
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
    , callPeakOptsMode  = Model
    , callPeakOptsCallSummits = False
    --, callPeakOptsBroad = False
    --, callPeakOptsBroadCutoff = 0.05
    }

-- | Call peaks using MACS2.
callPeaks :: SingI tags
          => FilePath                   -- ^ Ouptut file
          -> File tags 'Bed             -- ^ Sample
          -> Maybe (File tags 'Bed)    -- ^ Input/control sample
          -> CallPeakOptSetter          -- ^ Options
          -> IO (File '[] 'NarrowPeak)
callPeaks output target input setter = do
    macs2 output (target^.location) (fmap (^.location) input)
        fileFormat opt
    f <- frip (target^.location) output
    return $ location .~ output $ info .~ [("FRiP", T.pack $ show f)] $ emptyFile
  where
    opt = execState setter defaultCallPeakOpts
    fileFormat | target `hasTag` Pairend = "BEDPE"
               | otherwise = "BED"
{-# INLINE callPeaks #-}

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
            , "-g", T.pack $ opt^.gSize
            , "--outdir", T.pack tmp
            , "--tempdir", T.pack tmp
            , "--keep-dup", "all"
            , "-t", T.pack target ]
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
        mv (fromText $ T.pack $ tmp ++ "/NA_peaks.narrowPeak") $ fromText $
            T.pack output
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
        shelly $ mv (fromText $ T.pack final) $ fromText $ T.pack output
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
