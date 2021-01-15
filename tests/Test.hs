{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedStrings #-}
module Test
    ( tests
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline.NGS.Utils
import Bio.Pipeline.Barcode
import           Control.Lens
import           System.IO.Temp
import           Test.Tasty
import           Test.Tasty.Golden
import           Test.Tasty.HUnit

tests :: TestTree
tests = testGroup "Test"
    [ umiTest
    ]

filterBamTest :: TestTree
filterBamTest = goldenVsFile "Filter BAM Test" golden output $
    withSystemTempFile "tmp." $ \tmp h -> do
        filterBam "./" tmp (location .~ input $ emptyFile :: File '[PairedEnd] 'Bam) >>=
            sortBam "./" output
        return ()
  where
    golden = "tests/data/pairedend_filt.bam"
    output = "tests/data/pairedend_filt_test.bam"
    input = "tests/data/pairedend.bam"

umiTest :: TestTree
umiTest = testGroup "Barcode"
    [ testCase "conversion" $ map (intToDna . dnaToInt) barcodes @=? barcodes
    , testCase "mismatch" $ map (\(a,b,_) -> nMismatch (dnaToInt a) (dnaToInt b)) mis @=? map (\(_,_,x) -> x) mis
    ]
  where
    barcodes = ["AAAACCG", "ATCGGA", "CCCAATCC"]
    mis = [ ("AAAACCG", "AAAACCT", 1)
          , ("TAAACCG", "AAAACCT", 2)
          , ("TAAGCCG", "AAAACCT", 3)
          , ("AAAACCG", "TAAACCG", 1)
          , ("AAAACCG", "TTTTAAA", 7)
          , ("CCTGGAA", "CCTGGAA", 0)
          ]