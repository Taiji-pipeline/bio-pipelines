{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Bio.Pipeline.Demultiplex
    ( deBarcode
    , barcodeStat
    , dnaToInt
    , intToDna
    , genErrorBarcode1
    , nMismatch
    , mkBarcodeMap
    ) where

import Bio.Data.Fastq
import Conduit
import Data.List
import Data.List.Ordered
import Data.Ord
import qualified Data.ByteString                      as BS
import qualified Data.ByteString.Char8                as B
import qualified Data.IntMap.Strict as I
import qualified Data.IntSet as S
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
intToDna :: Int -> B.ByteString
intToDna = B.pack . reverse . go []
  where
    go !acc !x
        | m == 0 = c : acc
        | otherwise = go (c:acc) m
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

-- | Generate error barcodes within 1 edit distance to the input barcode.
genErrorBarcode1 :: Barcode -> [Barcode]
genErrorBarcode1 bc = map (\(n, i) -> bc + i * 5^n) $ filter ((/=0) . snd) $ go 0 bc
  where
    go !acc !x | x' == 0 = res
               | otherwise = res <> go (acc+1) x'
      where
        res = map (\i -> (acc, i - x `mod` 5)) [0..4]
        x' = x `div` 5
{-# INLINE genErrorBarcode1 #-}

nMismatch :: Barcode -> Barcode -> Int
nMismatch a b = go 0 a b
  where
    go !acc !x !y | m1 == 0 && m2 == 0 = acc'
                  | otherwise = go acc' m1 m2
      where
        acc' | x `mod` 5 == y `mod` 5 = acc
             | otherwise = acc + 1
        m1 = x `div` 5
        m2 = y `div` 5
{-# INLINE nMismatch #-}

mkBarcodeMap :: [Barcode] -> I.IntMap Barcode
mkBarcodeMap = I.fromList . concatMap f
  where
    f x = zip (x : genErrorBarcode1 x) $ repeat x
{-# INLINE mkBarcodeMap #-}

mkBarcodeMap' :: Double -> U.Vector (Barcode, Int) -> [([Int], Int)]
mkBarcodeMap' thres input = mkMap $ loop (I.fromList $ zip init $ repeat [])
    (S.fromList init) (S.fromList $ drop (truncate thres) $ map fst $ U.toList input)
  where
    mkMap m = map f $ filter (not . (`S.member` s)) init
      where
        s = S.fromList $ concatMap (m I.!) init
        f i = (nubSort $ go $ m I.! i, i)
          where
            go [] = []
            go xs = xs <> go (nubSort $ concatMap (\k -> I.findWithDefault [] k m) xs)
    countTable = I.fromList $ U.toList input
    init = take (truncate thres) $ map fst $ U.toList input
    loop acc whitelist candidate
        | S.null whitelist' = acc'
        | otherwise = loop acc' whitelist' candidate'
      where
        acc' = foldl' (\m (k,v) -> I.alter (f v) k m) acc assoc
          where
            f v Nothing = Just [v]
            f v (Just x) = Just $ v : x
        (assoc, whitelist', candidate') = findMatch whitelist candidate
    findMatch whitelist candidate = (assoc, whitelist', candidate')
      where
        assoc = concatMap (\x -> zip (filter (isConnected x) $ S.toList whitelist) $ repeat x) $ S.toList candidate
        whitelist' = S.fromList (map snd assoc) S.\\ whitelist
        candidate' = candidate S.\\ whitelist
    isConnected child parent = nMismatch child parent <= 1 && 
        countTable I.! parent >= 2 * countTable I.! child
