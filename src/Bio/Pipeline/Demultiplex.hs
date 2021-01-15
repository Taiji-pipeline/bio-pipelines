{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.Pipeline.Demultiplex
    ( deBarcode
    , barcodeStat
    , dnaToInt
    , intToDna
    --, mkBarcodeMap
    ) where

import Bio.Data.Fastq
import Conduit
import Data.Ord
import qualified Data.ByteString                      as BS
import qualified Data.ByteString.Char8                as B
import qualified Data.IntMap.Strict as I
import qualified Data.Vector.Algorithms.Intro as I
import qualified Data.Vector.Unboxed as U

type Barcode = Int

-- | Remove barcodes from the fastq record.
deBarcode :: Int   -- ^ Lengths of the barcodes
          -> Fastq
          -> Maybe (Barcode, Fastq)
deBarcode l Fastq{..}
    | nBelowQ 10 (B.take l fastqSeqQual) > 1 = Nothing
    | otherwise = Just (bc, fastq)
  where
    fastq = Fastq fastqSeqId (B.drop l fastqSeq) $ B.drop l fastqSeqQual
    bc = dnaToInt $ B.take l fastqSeq
{-# INLINE deBarcode #-}

-- | The total number of bases that below the quality threshold.
nBelowQ :: Int   -- ^ quality threshold
        -> BS.ByteString
        -> Int
nBelowQ q = BS.foldl' f 0
  where
    f acc x = if fromIntegral x - 33 < q then acc + 1 else acc
{-# INLINE nBelowQ #-}

barcodeStat :: MonadIO m
            => Int   -- ^ Length of the barcodes
            -> ConduitT Fastq o m (U.Vector (Int, Int))
barcodeStat l = U.modify (I.sortBy $ flip $ comparing snd) . U.fromList .
    I.toList <$> foldlC f I.empty
  where
    f m fastq = case deBarcode l fastq of
        Nothing -> m
        Just (bc, _) -> I.insertWith (+) bc 1 m

-- | Convert quinary DNA sequence (N: 0, A: 1, C: 2, G: 3, T: 4) to decimal.
dnaToInt :: B.ByteString -> Int
dnaToInt = fst . B.foldl' f (0, 1)
  where
    f (acc, i) 'N' = (acc, i * 5)
    f (acc, i) 'A' = (acc + 1 * i, i * 5)
    f (acc, i) 'C' = (acc + 2 * i, i * 5)
    f (acc, i) 'G' = (acc + 3 * i, i * 5)
    f (acc, i) 'T' = (acc + 4 * i, i * 5)
    f _ _          = error "Unexpected character!"
{-# INLINE dnaToInt #-}

-- | Convert decimal to quinary DNA sequence.
intToDna :: Int   -- ^ length of the resulting bytestring
         -> Int
         -> B.ByteString
intToDna n = B.pack . reverse . go 1 []
  where
    go !i !acc !x
        | m == 0 && i >= n = c : acc
        | otherwise = go (i+1) (c:acc) m
      where
        c = case x `mod` 5 of
            0 -> 'N'
            1 -> 'A'
            2 -> 'C'
            3 -> 'G'
            4 -> 'T'
            _ -> undefined
        m = x `div` 5
{-# INLINE intToDna #-}

{-
mkBarcodeMap :: M.HashMap BarCode Int -> M.HashMap BarCode BarCode
mkBarcodeMap = mkBarcodeMap' . mkUMINet

mkBarcodeMap' :: G.Gr (BarCode, Int) () -> M.HashMap BarCode BarCode
mkBarcodeMap' gr = M.fromList $
    concatMap (f . map (fromJust . G.lab gr)) $ G.components gr
  where
    f xs = let x = fst $ maximumBy (comparing snd) xs
           in zip (map fst xs) $ repeat x
{-# INLINE mkBarcodeMap' #-}

-- | Construct UMI network using directional approach.
mkUMINet :: M.HashMap BarCode Int -> G.Gr (BarCode, Int) ()
mkUMINet umi = G.mkGraph nds es
  where
    es = mapMaybe (uncurry isConnect) $ comb nds
    nds = zip [0..] $ M.toList umi
    isConnect (a, (x, c_x)) (b, (y, c_y))
        | nMismatch x y <= 1 = if c_x >= 2 * c_y - 1 
          then Just (a, b, ())
          else if c_y >= 2 * c_x - 1 then Just (b, a, ()) else Nothing
        | otherwise = Nothing
      where
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _ = []
{-# INLINE mkUMINet #-}

nMismatch :: BarCode -> BarCode -> Int
nMismatch a b = undefined

-}