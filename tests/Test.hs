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
    [ filterBamTest
    , umiTest
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
umiTest = testCase "UMI" $ sort res @=? sort ans
  where
    res = uniqBarcode umi
    ans = ["ACGT", "AAAT"]
    umi = [ ("ACGT", 456)
          , ("TCGT", 2)
          , ("AAAT", 90)
          , ("CCGT", 2)
          , ("ACAT", 72)
          , ("ACAG", 1) ]