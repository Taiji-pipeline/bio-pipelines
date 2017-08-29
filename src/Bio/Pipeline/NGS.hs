{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE TypeOperators     #-}

module Bio.Pipeline.NGS
    ( BWAOpts
    , BWAOptSetter
    , bwaCores
    , bwaSeedLen
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

    , mapFileWithDefName
    ) where

import           Bio.Data.Experiment
import           Control.Lens
import           Data.Bitraversable          (bitraverse)
import           Data.Promotion.Prelude.List (Delete, Elem, Insert)
import           Data.Singletons             (SingI)
import qualified Data.Text                   as T
import           Text.Printf                 (printf)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils

bwaAlign :: ( tags1' ~ Delete 'Gzip tags1
            , tags2' ~ Delete 'Gzip tags2
            , Elem 'Pairend tags2 ~ 'True )
         => (FilePath, String)
         -> FilePath
         -> BWAOptSetter
         -> Either (ATACSeq S (File tags1 'Fastq))
                   (ATACSeq S (File tags2 'Fastq, File tags2 'Fastq))
         -> IO ( Either (ATACSeq S (File tags1' 'Bam))
                        (ATACSeq S (File tags2' 'Bam)) )
bwaAlign (dir, suffix) idx opt = bitraverse fun1 fun2
  where
    fun1 = mapFileWithDefName (dir++"/") suffix $ \output ->
        bwaAlign_ output idx opt . Left
    fun2 = mapFileWithDefName (dir++"/") suffix $ \output ->
        bwaAlign_ output idx opt . Right

filterBam :: ( SingI tags, Experiment experiment
             , tags' ~ (Insert 'Sorted tags) )
          => (FilePath, String)
          -> experiment S (File tags 'Bam)
          -> IO (experiment S (File tags' 'Bam))
filterBam (dir, suffix) = mapFileWithDefName (dir++"/") suffix filterBam_

removeDuplicates :: (Experiment experiment, SingI tags)
                 => FilePath   -- ^ picard
                 -> (FilePath, String)
                 -> experiment S (File tags 'Bam)
                 -> IO (experiment S (File tags 'Bam))
removeDuplicates picard (dir, suffix) = mapFileWithDefName (dir++"/") suffix
    (removeDuplicates_ picard)

bamToBed :: ( Experiment experiment, SingI tags
            , Elem 'Sorted tags ~ 'True )
         => (FilePath, String)
         -> experiment S (File tags 'Bam)
         -> IO (experiment S (File (Insert 'Gzip tags) 'Bed))
bamToBed (dir, suffix) = mapFileWithDefName (dir ++ "/") suffix $ \output fl ->
    if fl `hasTag` Pairend
        then bam2Bed_ output (const True) fl
        else bam2BedPE_ output (const True) fl

concatBed :: (SingI tags, Experiment experiment)
          => (FilePath, String)
          -> experiment N (File tags 'Bed)
          -> IO (experiment S (File (Insert 'Gzip tags) 'Bed))
concatBed (dir, suffix) e = do
    fl <- concatBed_ output fls
    return $ e & replicates .~ return (Replicate fl [] 0)
  where
    fls = e^..replicates.folded.files
    output = printf "%s/%s_rep0%s" dir (T.unpack $ e^.eid) suffix

starAlign :: ( SingI tags1, SingI tags2, Elem 'Pairend tags2 ~ 'True
             , tags1' ~ Delete 'Gzip tags1
             , tags2' ~ Delete 'Gzip tags2 )
          => (FilePath, String)
          -> FilePath            -- ^ STAR genome index
          -> STAROptSetter       -- ^ Options
          -> Either (RNASeq S (File tags1 'Fastq))
                    (RNASeq S (File tags2 'Fastq, File tags2 'Fastq))
          -> IO ( Either (RNASeq S (File tags1' 'Bam, File tags1' 'Bam))
                         (RNASeq S (File tags2' 'Bam, File tags2' 'Bam)) )
starAlign (dir, suffix) idx setter = bitraverse fun1 fun2
  where
    fun1 = mapFileWithDefName (dir++"/") "" $ \output -> starAlign_
        (output ++ "_genome" ++ suffix) (output ++ "_anno" ++ suffix)
        idx setter . Left
    fun2 = mapFileWithDefName (dir++"/") "" $ \output -> starAlign_
        (output ++ "_genome" ++ suffix) (output ++ "_anno" ++ suffix)
        idx setter . Right

rsemQuant :: SingI tags
          => FilePath
          -> FilePath
          -> RSEMOptSetter
          -> RNASeq S (File tags 'Bam)
          -> IO (RNASeq S (File '[GeneQuant] 'Tsv, File '[TranscriptQuant] 'Tsv))
rsemQuant dir idx setter = mapFileWithDefName (dir++"/") "" $ \output ->
    rsemQuant_ output idx setter

mapFileWithDefName :: (Experiment experiment, Monad m, Traversable t)
                   => String   -- ^ Prefix
                   -> String   -- ^ Suffix
                   -> (FilePath -> file -> m file')
                   -> experiment t file
                   -> m (experiment t file')
mapFileWithDefName prefix suffix fn e = e & replicates.traverse %%~
    (\r -> r & files %%~ f r)
  where
    f r fl = fn output fl
      where
        output = printf "%s%s_rep%d%s" prefix (T.unpack $ e^.eid)
            (r^.number) suffix
{-# INLINE mapFileWithDefName #-}
