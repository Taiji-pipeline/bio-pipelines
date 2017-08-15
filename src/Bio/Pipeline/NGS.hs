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
    , bwaMaxMis
    , bwaReadTrim
    , defaultBWAOpts
    , bwaMkIndex
    , bwaAlign1
    , bwaAlign2
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

    , nameWith
    ) where

import           Bio.Data.Experiment
import           Control.Lens
import           Data.Either                 (fromLeft, fromRight)
import           Data.Promotion.Prelude.List (Delete, Elem, Insert)
import           Data.Singletons             (SingI)
import qualified Data.Text                   as T
import           Text.Printf                 (printf)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils

bwaAlign1 :: tags' ~ Delete 'Gzip tags
          => (FilePath, String)
          -> FilePath
          -> BWAOptSetter
          -> ATACSeq (File tags 'Fastq)
          -> IO (ATACSeq (File tags' 'Bam))
bwaAlign1 (dir, suffix) idx opt = nameWith (dir++"/") suffix $ \output input ->
    bwaAlign1_ output idx opt input

bwaAlign2 :: (tags' ~ Delete 'Gzip tags, Elem 'Pairend tags ~ 'True)
          => (FilePath, String)
          -> FilePath
          -> BWAOptSetter
          -> ATACSeq (File tags 'Fastq, File tags 'Fastq)
          -> IO (ATACSeq (File tags' 'Bam))
bwaAlign2 (dir, suffix) idx opt = nameWith (dir++"/") suffix $ \output (fw, rv) ->
    bwaAlign2_ output idx opt fw rv

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
                 -> IO (experiment (File tags 'Bam))
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
          -> IO (experiment (File (Insert 'Gzip tags) 'Bed))
concatBed (dir, suffix) e = do
    fl <- concatBed_ output fls
    return $ e & replicates .~ [ Replicate fl [] 0 ]
  where
    fls = e^..replicates.folded.files
    output = printf "%s/%s_rep0.%s" dir (T.unpack $ e^.eid) suffix

starAlign :: ( SingI tags1, SingI tags2, Elem 'Pairend tags2 ~ 'True
             , tags1' ~ Delete 'Gzip tags1
             , tags2' ~ Delete 'Gzip tags2 )
          => (FilePath, String)
          -> FilePath            -- ^ STAR genome index
          -> STAROptSetter       -- ^ Options
          -> Either (RNASeq (File tags1 'Fastq))
                    (RNASeq (File tags2 'Fastq, File tags2 'Fastq))
          -> IO ( Either (RNASeq (File tags1' 'Bam, File tags1' 'Bam))
                         (RNASeq (File tags2' 'Bam, File tags2' 'Bam))
                )
starAlign (dir, suffix) idx setter experiment = case experiment of
    Left e -> fmap Left $ e & replicates.traverse %%~ ( \r -> r & files %%~
        (fmap (fromLeft undefined) . fun (e^.eid) (r^.number) . Left) )
    Right e -> fmap Right $ e & replicates.traverse %%~ ( \r -> r & files %%~
        (fmap (fromRight undefined) . fun (e^.eid) (r^.number) . Right) )
  where
    fun id' rep fl = starAlign_ outputGenome outputAnno idx setter fl
      where
        outputGenome = printf "%s/%s_rep%d_genome.%s" dir (T.unpack id')
            rep suffix
        outputAnno = printf "%s/%s_rep%d_anno.%s" dir (T.unpack id')
            rep suffix

rsemQuant :: SingI tags
          => FilePath
          -> FilePath
          -> RSEMOptSetter
          -> RNASeq (File tags 'Bam)
          -> IO (RNASeq (File tags 'Tsv, File tags 'Tsv))
rsemQuant dir idx setter e = e & replicates.traverse %%~
    ( \r -> r & files %%~ fun (e^.eid) (r^.number) )
  where
    fun id' rep = rsemQuant_ outputPrefix idx setter
      where
        outputPrefix = printf "%s/%s_rep%d" dir (T.unpack id') rep

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
