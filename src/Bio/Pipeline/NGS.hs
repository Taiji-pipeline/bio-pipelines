{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE OverloadedLists      #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE TemplateHaskell      #-}
{-# LANGUAGE TypeOperators        #-}

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
    , concatBed

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

import           Bio.Data.Experiment
import           Control.Lens
import           Data.Promotion.Prelude.List (Delete, Elem, Insert)
import           Data.Singletons             (SingI)
import qualified Data.Text                   as T
import           Text.Printf                 (printf)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils

bwaAlign :: (tags1' ~ Delete 'Gzip tags1, tags2' ~ Delete 'Gzip tags2)
         => (FilePath, String)
         -> FilePath
         -> BWAOptSetter
         -> ATACSeq (MaybePair tags1 tags2 'Fastq)
         -> IO (ATACSeq (EitherTag tags1' tags2' 'Bam))
bwaAlign (dir, suffix) idx opt = nameWith (dir++"/") suffix $ \output input ->
    bwaAlign_ output idx opt input

filterBam :: ( SingI tags, Experiment experiment
             , tags' ~ (Insert 'Sorted tags) )
          => (FilePath, String)
          -> experiment (File tags 'Bam)
          -> IO (experiment (File tags' 'Bam))
filterBam (dir, suffix) = nameWith (dir++"/") suffix filterBam_

removeDuplicates :: (Experiment experiment, SingI tags)
                 => FilePath   -- ^ picard
                 -> (FilePath, String)
                 -> experiment (File tags 'Bam)
                 -> IO (experiment (File tags 'Bam, File '[] 'Other))
removeDuplicates picard (dir, suffix) = nameWith (dir++"/") suffix
    (removeDuplicates_ picard)

bamToBed :: ( Experiment experiment, SingI tags
            , Elem 'Sorted tags ~ 'True )
         => (FilePath, String)
         -> experiment (File tags 'Bam)
         -> IO (experiment (File (Insert 'Gzip tags) 'Bed))
bamToBed (dir, suffix) = nameWith (dir ++ "/") suffix fn
  where
    fn output fl = if fl `hasTag` Pairend
        then bam2Bed_ output (const True) fl
        else bam2BedPE_ output (const True) fl

concatBed :: (SingI tags, Experiment experiment)
          => (FilePath, String)
          -> experiment (File tags 'Bed)
          -> IO (experiment (File '[Gzip] 'Bed))
concatBed (dir, suffix) e = do
    fl <- concatBed_ output fls
    return $ e & replicates .~ [ Replicate fl [] 0 ]
  where
    fls = e^..replicates.folded.files
    output = printf "%s/%s_rep0.%s" dir (T.unpack $ e^.eid) suffix

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
