{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}
module Bio.Pipeline.Barcode
    ( deBarcode
    , barcodeStat
    ) where

import Bio.Data.Fastq
import Conduit
import Data.List
import qualified Data.ByteString                      as BS
import qualified Data.ByteString.Char8                as B
import qualified Data.HashMap.Strict as M

barcodeStat :: Monad m
            => [Int]   -- ^ Lengths of the barcodes
            -> ConduitT Fastq o m [M.HashMap Int Int]
barcodeStat lengths = foldlC f $ replicate (length lengths) M.empty
  where
    f m fastq = case deBarcode lengths fastq of
        Nothing -> m
        Just (barcodes, _) -> zipWith (\b -> M.insertWith (+) b 1) barcodes m

deBarcode :: [Int]   -- ^ Lengths of the barcodes
          -> Fastq
          -> Maybe ([Int], Fastq)
deBarcode lengths Fastq{..}
    | any (\x -> nBelowQ 10 x > 1) (go lengths fastqSeqQual) = Nothing
    | otherwise = Just (barcodes, fastq)
  where
    fastq = Fastq ident (B.drop n fastqSeq) $ B.drop n fastqSeqQual
    ident = B.intercalate ":" (map (B.pack . show) barcodes) <>
        "_barcode_" <> fastqSeqId
    barcodes = map dnaToInt $ go lengths fastqSeq
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
