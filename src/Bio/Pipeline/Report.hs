{-# LANGUAGE DataKinds #-}
{-# LANGUAGE BangPatterns #-}
module Bio.Pipeline.Report where

import Conduit
import Control.Lens ((^.))
import Bio.HTS
import Bio.Data.Experiment

bamStat :: File tags 'Bam -> IO (Int, Int)
bamStat input = withBamFile (input^.location) $ \h ->
    runConduit $ readBam h .| foldlC f (0,0)
  where
    f (!mapped, !total) bam =
        (if isUnmapped (flag bam) then mapped else mapped + 1, total + 1)
