{-# LANGUAGE DataKinds #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Bio.Pipeline.Report where

import Conduit
import qualified Data.ByteString.Char8 as B
import Control.Lens ((^.))
import qualified Data.Map.Strict as M
import Bio.HTS
import Bio.Data.Experiment

data BamStats = BamStats
    { _mapped_tags :: Int
    , _total_tags :: Int
    , _mapped_tags_by_chrom :: [(B.ByteString, Int)]
    }

bamStat :: File tags 'Bam -> IO BamStats
bamStat input = do
    header <- getBamHeader $ input^.location
    (mapped, total) <- runResourceT $ runConduit $
        streamBam (input^.location) .| foldMC (f header) (M.empty,0)
    return $ BamStats (sum $ M.elems mapped) total $ M.toList mapped
  where
    f header (!mapped, !total) bam
        | isUnmapped (flag bam) = return (mapped, total + 1)
        | otherwise = return $ case refName header bam of
            Nothing -> (mapped, total + 1)
            Just chr -> (M.alter g chr mapped, total + 1)
    g Nothing = Just 1
    g (Just x) = Just $ x + 1
