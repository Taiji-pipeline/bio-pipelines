module Bio.Pipeline 
    ( module Bio.Pipeline.Download
    , module Bio.Pipeline.CallPeaks
    , module Bio.Pipeline.Barcode
    , module Bio.Pipeline.Report
    , module Bio.Pipeline.Utils
    , module Bio.Pipeline.NGS.BWA
    , module Bio.Pipeline.NGS.RSEM
    , module Bio.Pipeline.NGS.STAR
    , module Bio.Pipeline.NGS.Utils
    , RAWInput
    , FASTQInput
    , getFastq
    ) where

import Bio.Data.Experiment

import Bio.Pipeline.Download
import Bio.Pipeline.CallPeaks
import Bio.Pipeline.Barcode
import Bio.Pipeline.Report
import Bio.Pipeline.Utils
import Bio.Pipeline.NGS.BWA
import Bio.Pipeline.NGS.RSEM
import Bio.Pipeline.NGS.STAR
import Bio.Pipeline.NGS.Utils

type RAWInput e = e N [Either SomeFile (SomeFile, SomeFile)]

type FASTQInput e = [e S (Either (SomeTags 'Fastq) (SomeTags 'Fastq, SomeTags 'Fastq))]

getFastq :: Experiment e
         => [RAWInput e]
         -> [FASTQInput e]
getFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap castFile (bimap castFile castFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq