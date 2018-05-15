{-# LANGUAGE DataKinds #-}
module Test
    ( tests
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline.NGS.Utils
import           Control.Lens
import           System.IO.Temp
import           Test.Tasty
import           Test.Tasty.Golden

tests :: TestTree
tests = testGroup "Test"
    [ filterBamTest
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
