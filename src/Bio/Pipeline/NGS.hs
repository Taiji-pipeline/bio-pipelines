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
    , removeDuplicates
    , bamToBed
    , concatBed

    , mapFileWithDefName
    ) where

import           Bio.Data.Experiment
import           Control.Lens
import           Data.Bitraversable          (bitraverse)
import           Data.Promotion.Prelude.List (Delete, Elem)
import           Data.Singletons             (SingI)
import qualified Data.Text                   as T
import           Text.Printf                 (printf)

import           Bio.Pipeline.NGS.BWA
import           Bio.Pipeline.NGS.Utils

bwaAlign :: ( tags1' ~ Delete 'Gzip tags1
            , tags2' ~ Delete 'Gzip tags2
            , Elem 'PairedEnd tags2 ~ 'True )
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

removeDuplicates :: (Experiment experiment, SingI tags)
                 => FilePath   -- ^ picard
                 -> (FilePath, String)
                 -> experiment S (File tags 'Bam)
                 -> IO (experiment S (File tags 'Bam))
removeDuplicates picard (dir, suffix) = mapFileWithDefName (dir++"/") suffix
    (removeDuplicates_ picard)

bamToBed :: ( Experiment experiment, SingI tags
            , Elem 'CoordinateSorted tags ~ 'True )
         => (FilePath, String)
         -> experiment S (File tags 'Bam)
         -> IO (experiment S (File (Insert' 'Gzip tags) 'Bed))
bamToBed (dir, suffix) = mapFileWithDefName (dir ++ "/") suffix $ \output fl ->
    if fl `hasTag` PairedEnd
        then bam2Bed_ output (const True) fl
        else bam2BedPE_ output (const True) fl

concatBed :: ( Experiment experiment
             , Elem 'Gzip tags1 ~ 'False
             , Elem 'Gzip tags2 ~ 'True )
          => (FilePath, String)
          -> experiment N (Either (File tags1 'Bed) (File tags2 'Bed))
          -> IO (experiment S (File '[Gzip] 'Bed))
concatBed (dir, suffix) e = do
    fl <- concatBed_ output fls
    return $ e & replicates .~ return (Replicate fl [] 0)
  where
    fls = e^..replicates.folded.files
    output = printf "%s/%s_rep0%s" dir (T.unpack $ e^.eid) suffix

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
