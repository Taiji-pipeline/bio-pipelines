{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE ExtendedDefaultRules #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE OverloadedLists      #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE TemplateHaskell      #-}

module Bio.Pipeline.NGS
    ( BWAOpts
    , BWAOptSetter
    , bwaCores
    , bwaSeedLen
    , bwaMaxMis
    , bwaReadTrim
    , defaultBWAOpts
    , bwaMkIndex
    , bwaAlign
    , filterBam
    --, removeDuplicates
    --, bam2Bed
    --, bam2BedPE
    , mergeReplicatesBed

    , STAROpts
    , STAROptSetter
    , starCmd
    , starCores
    , starSort
    , starTmpDir
    , starMkIndex
    , starAlign

    , RSEMOpts
    , RSEMOptSetter
    , rsemPath
    , rsemCores
    , rsemSeed
    , rsemMkIndex
    , rsemQuant
    ) where

import           Bio.Data.Bam             (bamToBed, readBam, runBam,
                                           sortedBamToBedPE)
import           Bio.Data.Bed             (BED, BED3 (..), BEDLike (..), toLine)
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib        (gzip, ungzip)
import           Data.Maybe               (fromJust)
import Data.Typeable (Typeable)
import           Data.Tagged              (Tagged (..), untag)
import qualified Data.Text                as T
import           Shelly                   (escaping, fromText, mv, run_, shelly,
                                           silently)
import           System.IO.Temp           (withTempDirectory)
import Text.Printf (printf)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils


bwaAlign :: FilePath
         -> FilePath
         -> BWAOptSetter
         -> ATACSeq (MaybeTagged 'GZipped (MaybePaired (File 'Fastq)))
         -> IO (ATACSeq (MaybeTagged 'Pairend (File 'Bam)))
bwaAlign dir index opt atac = atac & replicates.traverse %%~
    (\r -> r & files %%~ align r)
  where
    align r = bwaAlign_ output index opt . fst . untagMaybe
      where
        output = printf "%s/%s_rep%d.bam" dir (T.unpack $ atac^.eid) (r^.number)

filterBam :: MayHave 'Pairend tags
          => FilePath
          -> ATACSeq (Tagged tags (File 'Bam))
          -> IO (ATACSeq (Tagged tags (File 'Bam)))
filterBam dir e = e & replicates.traverse %%~
    (\r -> r & files %%~ fn r)
  where
    fn r fl = filterBam_ output fl
      where
        output = printf "%s/%s_rep%d.filt.bam" dir (T.unpack $ e^.eid) (r^.number)

        {-
removeDuplicates :: FilePath
                 -> ATACSeq (MaybeTagged Pairend (Tagged Sorted (File 'Bam)))
                 -> IO (ATACSeq (MaybeTagged Pairend (File 'Bam)))
                 -}
