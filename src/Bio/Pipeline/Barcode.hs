{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.Pipeline.Barcode
    ( BarCode(..)
    , deBarcode
    , barcodeStat
    , mkBarcodeMap
    ) where

import Bio.Data.Fastq
import Conduit
import Data.List
import Data.Ord
import qualified Data.ByteString                      as BS
import qualified Data.ByteString.Char8                as B
import qualified Data.HashMap.Strict as M
import qualified Data.Graph.Inductive as G
import Data.String (IsString)
import Data.Maybe
import Data.Hashable (Hashable)

newtype BarCode = BarCode { unBarCode :: B.ByteString }
      deriving (Eq, Ord, Hashable, IsString, Show)

barcodeStat :: Monad m
            => [Int]   -- ^ Lengths of the barcodes
            -> ConduitT Fastq o m [M.HashMap BarCode Int]
barcodeStat lengths = foldlC f $ replicate (length lengths) M.empty
  where
    f m fastq = case deBarcode lengths fastq of
        Nothing -> m
        Just (barcodes, _) -> zipWith (\b -> M.insertWith (+) b 1) barcodes m

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
        | go 0 0 = if c_x >= 2 * c_y - 1 
          then Just (a, b, ())
          else if c_y >= 2 * c_x - 1 then Just (b, a, ()) else Nothing
        | otherwise = Nothing
      where
        n = B.length $ unBarCode x
        go :: Int -> Int -> Bool
        go !mis !i
            | mis > 1 = False
            | i >= n = True
            | otherwise = if unBarCode x `B.index` i == unBarCode y `B.index` i
              then go mis $ i+1 
              else go (mis+1) $ i+1
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _ = []
{-# INLINE mkUMINet #-}

-- | Remove barcodes from the fastq record.
deBarcode :: [Int]   -- ^ Lengths of the barcodes
          -> Fastq
          -> Maybe ([BarCode], Fastq)
deBarcode lengths Fastq{..}
    | any (\x -> nBelowQ 10 x > 1) (go lengths fastqSeqQual) = Nothing
    | otherwise = Just (map BarCode barcodes, fastq)
  where
    fastq = Fastq fastqSeqId (B.drop n fastqSeq) $ B.drop n fastqSeqQual
    barcodes = go lengths fastqSeq
    go (x:xs) s = B.take x s : go xs (B.drop x s)
    go _ _ = []
    n = foldl1' (+) lengths
{-# INLINE deBarcode #-}

-- | The total number of bases that below the quality threshold.
nBelowQ :: Int   -- ^ quality threshold
        -> BS.ByteString
        -> Int
nBelowQ q = BS.foldl' f 0
  where
    f acc x = if fromIntegral x - 33 < q then acc + 1 else acc
{-# INLINE nBelowQ #-}

{-
-- | Convert quinary DNA sequence (A: 0, C: 1, G: 2, T: 3, N: 4) to decimal.
dnaToInt :: B.ByteString -> Int
dnaToInt = fst . B.foldl' f (0, 1)
  where
    f (acc, i) 'A' = (acc, i * 5)
    f (acc, i) 'C' = (acc + 1 * i, i * 5)
    f (acc, i) 'G' = (acc + 2 * i, i * 5)
    f (acc, i) 'T' = (acc + 3 * i, i * 5)
    f (acc, i) 'N' = (acc + 4 * i, i * 5)
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
            0 -> 'A'
            1 -> 'C'
            2 -> 'G'
            3 -> 'T'
            4 -> 'N'
            _ -> undefined
        m = x `div` 5
{-# INLINE intToDna #-}
-}
