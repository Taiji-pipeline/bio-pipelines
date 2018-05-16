{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE TypeOperators     #-}

module Bio.Pipeline.NGS
    ( bwaMkIndex
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
