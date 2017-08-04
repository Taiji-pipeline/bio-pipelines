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
    , removeDuplicates
    , bamToBed

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

import           Bio.Data.Bed             (BED, BED3 (..), BEDLike (..), toLine)
import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib        (gzip, ungzip)
import           Data.Maybe               (fromJust)
import Data.Typeable (Typeable)
import qualified Data.Text                as T
import           Shelly                   (escaping, fromText, mv, run_, shelly,
                                           silently)
import Text.Printf (printf)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils

bwaAlign :: (MayHave 'Pairend tags, MayHave 'GZipped tags)
         => (FilePath, String)
         -> FilePath
         -> BWAOptSetter
         -> ATACSeq (MaybePaired (File tags 'Fastq))
         -> IO (ATACSeq (File (Remove 'GZipped tags) 'Bam))
bwaAlign (dir, suffix) index opt = nameWith (dir++"/") suffix $ \output input ->
    bwaAlign_ output index opt input

filterBam :: ( MayHave 'Pairend tags, Experiment experiment
             , tags' ~ Insert 'Sorted tags )
          => (FilePath, String)
          -> experiment (File tags 'Bam)
          -> IO (experiment (File tags' 'Bam))
filterBam (dir, suffix) = nameWith (dir++"/") suffix filterBam_

removeDuplicates :: (Experiment experiment, MayHave 'Pairend tags)
                 => FilePath   -- ^ picard
                 -> (FilePath, String)
                 -> experiment (File tags 'Bam)
                 -> IO (experiment (File tags 'Bam, File '[] 'Other))
removeDuplicates picard (dir, suffix) = nameWith (dir++"/") suffix
    (removeDuplicates_ picard)

bamToBed :: ( Experiment experiment
            , MayHave 'Pairend tags, Elem 'Sorted tags ~ 'True)
         => (FilePath, String)
         -> experiment (File tags 'Bam)
         -> IO (experiment (File (Insert 'GZipped tags) 'Bed))
bamToBed (dir, suffix) = nameWith (dir ++ "/") suffix fn
  where
    fn output fl = if isPairend fl
        then bam2Bed_ output (const True) fl
        else bam2BedPE_ output (const True) fl

nameWith :: (Experiment experiment, Monad m)
         => String   -- ^ Prefix
         -> String   -- ^ Suffix
         -> (FilePath -> file -> m file')
         -> experiment file
         -> m (experiment file')
nameWith prefix suffix fn e = e & replicates.traverse %%~ (\r -> r & files %%~ f r)
  where
    f r fl = fn output fl
      where
        output = printf "%s%s_rep%d.%s" prefix (T.unpack $ e^.eid)
            (r^.number) suffix
{-# INLINE nameWith #-}
