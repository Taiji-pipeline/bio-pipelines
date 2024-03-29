{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE OverloadedLists      #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE ScopedTypeVariables  #-}
{-# LANGUAGE TemplateHaskell      #-}
{-# LANGUAGE TypeFamilies         #-}
{-# LANGUAGE TypeOperators        #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
module Bio.Pipeline.NGS.Utils
    ( filterBam
    , filterBamSort
    , sortBam
    , sortBamByName
    , removeDuplicates
    , bam2Bed
    , bam2BedPE
    , concatBed
    , bedToBigBed
    , bedToBigWig
    , bedToBigWigC
    , bedGraphToBigWig
    , bamToUniqFragment
    ) where

import           Bio.Data.Bam                (bamToBedC, sortedBamToBedPE)
import           Bio.Data.Bed
import           Bio.Data.Bed.Utils (clipBed)
import           Bio.Data.Experiment
import           Conduit
import Data.Conduit.Internal (zipSinks)
import qualified Data.Conduit.List as CL
import           Control.Lens
import           Control.Monad.State.Lazy
import           Data.Conduit.Zlib           (gzip, ungzip)
import           Data.Maybe                  (fromJust)
import           Prelude.Singletons (Elem, If, SingI)
import Data.List.Singletons (Delete)
import qualified Data.ByteString.Char8 as B
import qualified Data.Text                   as T
import           Shelly  hiding (FilePath)
import qualified Data.HashMap.Strict as M
import qualified Data.Map.Strict as Map
import Data.Double.Conversion.ByteString (toShortest)
import Bio.HTS.Utils (makeKeyPair, BAMKey(..))
import Bio.HTS.Types
import Bio.HTS.BAM
import Data.Char (isDigit)
import Data.Function (on)
import Data.ByteString.Lex.Integral (packDecimal)

import Bio.Pipeline.Utils

-- | Remove low quality and redundant tags, fill in mate information.
-- Use up to two threads.
filterBam :: ( SingI tags, tags' ~ If (Elem PairedEnd tags)
               (Insert' 'NameSorted (Delete 'CoordinateSorted tags)) tags )
          => FilePath  -- ^ temp dir
          -> FilePath  -- ^ output
          -> File tags 'Bam
          -> IO (File tags' 'Bam)
filterBam tmpDir output fl = withTempDir (Just tmpDir) $ \tmp -> do
    let input = T.pack $ fl^.location
        tmp_sort = T.pack $ tmp ++ "/tmp_sort"
    shelly $ escaping False $ silently $ if isPair
        then bashPipeFail bash_ "samtools"
            [ "view", "-@", "2", "-f", "2", "-F", "0x70c", "-q", "30", "-u", input, "|"
            , "samtools", "sort", "-@", "2", "-", "-n", "-T", tmp_sort, "-m", "4G", "-l", "0", "|" 
            , "samtools", "fixmate", "-@", "2", "-r", "-m", "-", "-", "|"
            , "samtools", "view", "-@", "2", "-F", "1804", "-f", "2", "-b", "-", ">"
            , T.pack output ]
        else run_ "samtools" [ "view", "-@", "2", "-F", "0x70c", "-q", "30", "-b", input
            , ">", T.pack output ]
    return $ location .~ output $ emptyFile
  where
    isPair = fl `hasTag` PairedEnd
{-# INLINE filterBam #-}

-- | Remove low quality and redundant tags, fill in mate information.
-- Finally, sort the bam by coordinate.
filterBamSort :: ( SingI tags
                 , tags' ~ Insert' 'CoordinateSorted (Delete 'NameSorted tags) )
              => FilePath  -- ^ temp dir
              -> FilePath  -- ^ output
              -> File tags 'Bam
              -> IO (File tags' 'Bam)
filterBamSort tmpDir output fl = withTempDir (Just tmpDir) $ \tmp -> do
    let input = escapeFileName $ T.pack $ fl^.location
    shelly $ silently $ if isPair
        then bashPipeFail bash_ "samtools"
            [ "view", "-f", "2", "-F", "0x70c", "-q", "30", "-u", input, "|"
            , "samtools", "sort", "-@", "2", "-", "-n", "-m", "4G", "-l", "0",
                "-T", escapeFileName $ T.pack tmp <> "/nsrt", "|"
            , "samtools", "fixmate", "-r", "-m", "-", "-", "|"
            , "samtools", "view", "-F", "1804", "-f", "2", "-u", "-", "|"
            , "samtools", "sort", "-@", "2", "-", "-T", T.pack tmp <> "/csrt",
                "-m", "4G", "-l", "9", "-o", escapeFileName $ T.pack output ]
        else bashPipeFail bash_ "samtools"
            [ "view", "-F", "0x70c", "-q", "30", "-u", input, "|"
            , "samtools", "sort", "-@", "2", "-", "-T", escapeFileName $ T.pack tmp <> "/csrt",
                "-m", "4G", "-l", "9", "-o", escapeFileName $ T.pack output ]
    return $ location .~ output $ emptyFile
  where
    isPair = fl `hasTag` PairedEnd
{-# INLINE filterBamSort #-}

sortBam :: ( SingI tags, tags' ~ Insert' 'CoordinateSorted
             (Delete 'NameSorted tags) )
        => FilePath    -- ^ temp dir
        -> FilePath    -- ^ output
        -> File tags 'Bam
        -> IO (File tags' 'Bam)
sortBam tmpDir output fl = withTempDir (Just tmpDir) $ \tmp -> do
    let input = T.pack $ fl^.location
        tmp_sort = T.pack $ tmp ++ "/tmp_sort"
    shelly $ silently $ run_ "samtools"
        [ "sort", "-@", "2", input, "-T", tmp_sort, "-l", "9", "-o", T.pack output ]
    return $ location .~ output $ emptyFile
{-# INLINE sortBam #-}

sortBamByName :: ( SingI tags, tags' ~ Insert' 'NameSorted
                   (Delete 'CoordinateSorted tags) )
              => FilePath    -- ^ temp dir
              -> FilePath    -- ^ output
              -> File tags 'Bam
              -> IO (File tags' 'Bam)
sortBamByName tmpDir output fl = withTempDir (Just tmpDir) $ \tmp -> do
    let input = T.pack $ fl^.location
        tmp_sort = T.pack $ tmp ++ "/tmp_sort"
    shelly $ silently $ run_ "samtools"
        [ "sort", "-@", "2", input, "-n", "-T", tmp_sort, "-l", "9", "-o", T.pack output ]
    return $ location .~ output $ emptyFile
{-# INLINE sortBamByName #-}

-- | Remove duplicates
removeDuplicates :: Elem 'CoordinateSorted tags ~ 'True
                 => FilePath
                 -> File tags 'Bam
                 -> IO (File tags 'Bam)
removeDuplicates output input = shelly $ do
    _ <- run "samtools" [ "markdup", T.pack $ input^.location
        , T.pack output, "-r", "-s" ]
    qc <- Map.fromList . parseQC <$> lastStderr
    return $ info .~ qc $ location .~ output $ emptyFile
  where
    parseQC = concatMap (parse . T.filter (/=':')) . filter (not . T.null) . T.lines
      where
        parse x =
            let (key, txt) = T.break isDigit x
                (val, rest) = T.span isDigit txt
            in if T.null val
                then []
                else (T.toUpper $ T.strip key, val) : parse rest
{-# INLINE removeDuplicates #-}

bam2Bed :: FilePath
        -> (BED -> Bool)  -- ^ Filtering function
        -> File tags 'Bam -> IO (File (Insert' 'Gzip tags) 'Bed)
bam2Bed output fn fl = do
    header <- getBamHeader $ fl^.location
    runResourceT $ runConduit $ streamBam (fl^.location) .| bamToBedC header .|
        filterC fn .| mapC toLine .| unlinesAsciiC .| gzip .| sinkFileBS output
    return $ location .~ output $ emptyFile
{-# INLINE bam2Bed #-}

-- | Convert name sorted BAM to BEDPE suitable for MACS2.
bam2BedPE :: Elem 'NameSorted tags ~ 'True
          => String
          -> ((BED, BED) -> Bool)
          -> File tags 'Bam
          -> IO (File (Insert' 'Gzip tags) 'Bed)
bam2BedPE output fn fl = do
    header <- getBamHeader $ fl^.location
    runResourceT $ runConduit $ streamBam (fl^.location) .|
        sortedBamToBedPE header .| filterC fn .| concatMapC f .| mapC toLine .|
        unlinesAsciiC .| gzip .| sinkFileBS output
    return $ location .~ output $ emptyFile
  where
    f (b1, b2)
        | b1^.chrom /= b2^.chrom || b1^.strand == b2^.strand = Nothing
        | otherwise =
            let left = if fromJust (b1^.strand)
                    then b1^.chromStart else b2^.chromStart
                right = if not (fromJust $ b2^.strand)
                    then b2^.chromEnd else b1^.chromEnd
            in if left < right
                  then Just (asBed (b1^.chrom) left right :: BED3)
                  else error "Left coordinate is larger than right coordinate."
{-# INLINE bam2BedPE #-}

-- | Merge multiple BED files.
concatBed :: (Elem 'Gzip tags1 ~ 'False, Elem 'Gzip tags2 ~ 'True)
          => FilePath
          -> [Either (File tags1 'Bed) (File tags2 'Bed)]
          -> IO (File '[Gzip] 'Bed)
concatBed output fls = do
    runResourceT $ runConduit $ source .| gzip .| sinkFile output
    return $ location .~ output $ emptyFile
  where
    source = forM_ fls $ \fl -> case fl of
        Left fl'  -> sourceFileBS (fl'^.location)
        Right fl' -> sourceFileBS (fl'^.location) .| ungzip
{-# INLINE concatBed #-}

-- | Create a bigbed file from a bed file.
bedToBigBed :: Elem 'Gzip tags ~ 'False
            => FilePath    -- ^ Output
            -> [(B.ByteString, Int)]   -- ^ Chromosome sizes
            -> File tags 'Bed
            -> IO (File tags 'BigBed)
bedToBigBed output chrSizes input = shelly $ test_px "bedToBigBed" >>= \case
    False -> error "Please download: bedToBigBed"
    True -> withTempDir (Just "./") $ \dir -> do
        let tmpSort = T.pack $ dir ++ "/tmp_sort.bed"
            tmpChr = T.pack $ dir ++ "/tmp_chr.txt"
        escaping False $ run_ "sort" ["-T", T.pack dir, "-k1,1", "-k2,2n"
            , T.pack $ input^.location, ">", tmpSort]
        liftIO $ B.writeFile (T.unpack tmpChr) $ B.unlines $
            map (\(a,b) -> a <> "\t" <> B.pack (show b)) chrSizes
        run_ "bedToBigBed" [tmpSort, tmpChr, T.pack output]
        return $ location .~ output $ emptyFile
{-# INLINE bedToBigBed #-}

-- | Create a bigwig file from a bedgraph file.
bedGraphToBigWig :: FilePath   -- output
                 -> [(B.ByteString, Int)]   -- ^ Chromosome sizes
                 -> [BED3]  -- ^ Blacklist regions
                 -> FilePath
                 -> FilePath  -- ^ Input
                 -> IO ()
bedGraphToBigWig output chrSizes blacklist tmpdir input = shelly (test_px "bedGraphToBigWig") >>= \case
    False -> error "Please download: bedGraphToBigWig"
    True -> withTempDir (Just tmpdir) $ \dir -> do
        let tmp1 = dir ++ "/tmp1"
            tmp2 = dir ++ "/tmp2"
            tmpChr = dir ++ "/chr"
        B.writeFile tmpChr $ B.unlines $
            map (\(a,b) -> a <> "\t" <> B.pack (show b)) chrSizes
        runResourceT $ runConduit $ streamBed input .|
            filterC (\x -> not $ isIntersected blacklist' (x :: BEDGraph)) .|
            clipBed chrSizes .| sinkFileBed tmp1
        shelly $ do
            setenv "LC_COLLATE" "C"
            escaping False $ run_ "sort"
                ["-S", "4G", "-k", "1,1", "-k2,2n", T.pack tmp1, ">", T.pack tmp2]
        shelly $ run_ "bedGraphToBigWig" [T.pack tmp2, T.pack tmpChr, T.pack output]
  where
    blacklist' = bedToTree const $ zip blacklist $ repeat ()
{-# INLINE bedGraphToBigWig #-}

-- | Create a bigwig file from a bed file.
bedToBigWig :: Elem 'Gzip tags ~ 'True
            => FilePath   -- output
            -> [(B.ByteString, Int)]   -- ^ Chromosome sizes
            -> [BED3]  -- ^ Blacklist regions
            -> FilePath
            -> File tags 'Bed
            -> IO ()
bedToBigWig output chrSizes blacklist tmpdir input = shelly (test_px "bedGraphToBigWig") >>= \case
    False -> error "Please download: bedGraphToBigWig"
    True -> withTempDir (Just tmpdir) $ \dir -> do
        let tmp1 = dir ++ "/tmp1"
            tmp2 = dir ++ "/tmp2"
            tmpChr = dir ++ "/chr"
        B.writeFile tmpChr $ B.unlines $
            map (\(a,b) -> a <> "\t" <> B.pack (show b)) chrSizes

        numReads <- extendBed tmp1 chrSizes $ input^.location
        shelly $ do
            setenv "LC_COLLATE" "C"
            escaping False $ run_ "sort"
                ["-S", "4G", "-k", "1,1", "-k2,2n", T.pack tmp1, ">", T.pack tmp2]
        mkBedGraph tmp1 tmp2 numReads
        shelly $ run_ "bedGraphToBigWig" [T.pack tmp1, T.pack tmpChr, T.pack output]
  where
    blacklist' = bedToTree const $ zip blacklist $ repeat ()
    extendBed out chr fl = do
        (n, _) <- runResourceT $ runConduit $ streamBedGzip fl .|
            filterC (not . isIntersected blacklist') .| concatMapC f .|
            zipSinks (mapC size .| sumC) (sinkFileBed out)
        return n
      where
        f :: BED -> Maybe BED3
        f bed = case M.lookup (bed^.chrom) chrSize of
            Nothing -> Nothing
            Just n -> Just $ case bed^.strand of
                Just False -> BED3 (bed^.chrom) (max 0 $ bed^.chromEnd - 100) (bed^.chromEnd)
                _ -> BED3 (bed^.chrom) (bed^.chromStart) (min n $ bed^.chromStart + 100)
        chrSize = M.fromList chr
{-# INLINE bedToBigWig #-}

-- | Create a bigwig file from a bed file.
bedToBigWigC :: BEDLike b
             => FilePath   -- output
             -> FilePath   -- tmp
             -> [(B.ByteString, Int)]   -- ^ Chromosome sizes
             -> ConduitT b Void (ResourceT IO) ()
bedToBigWigC output tmpdir chrSizes = do
    liftIO $ shelly (test_px "bedGraphToBigWig") >>= \case
        False -> error "Please download: bedGraphToBigWig"
        True -> return ()
    (numReads, tmpBed) <- zipSinks (mapC size .| sumC) $
        mapC ((toLine :: BED3 -> B.ByteString) . convert) .| unlinesAsciiC .| 
        sinkTempFile tmpdir "tmp.bed."
    tmpChr <- yieldMany chrSizes .| mapC (\(a,b) -> a <> "\t" <> B.pack (show b)) .|
        unlinesAsciiC .| sinkTempFile tmpdir "tmp.chr."
    tmpF <- yield "" .| sinkTempFile tmpdir "tmp.file."
    liftIO $ shelly $ do
        setenv "LC_COLLATE" "C"
        escaping False $ run_ "sort"
            ["-S", "4G", "-k", "1,1", "-k2,2n", T.pack tmpBed, ">", T.pack tmpF]
    liftIO $ do
        mkBedGraph tmpBed tmpF numReads
        shelly $ run_ "bedGraphToBigWig" [T.pack tmpBed, T.pack tmpChr, T.pack output]
{-# INLINE bedToBigWigC #-}

mkBedGraph :: FilePath  -- ^ Output
           -> FilePath  -- ^ Coordinate sorted bed files
           -> Int
           -> IO ()
mkBedGraph output input nReads = runResourceT $ runConduit $ streamBed input .|
    mergeSortedBedWith (countOverlapped :: [BED3] -> [(BED3, Int)]) .|
    concatC .| mapC f .| unlinesAsciiC .| sinkFile output
  where
    f (bed, x) = toLine bed <> "\t" <> toShortest (fromIntegral x / n)
    n = fromIntegral nReads / 1000000 / 1000
{-# INLINE mkBedGraph #-}

escapeFileName :: T.Text -> T.Text
escapeFileName = T.replace "'" "'\"\\'\"'"
{-# INLINE escapeFileName #-}

bamToUniqFragment :: FilePath
                  -> BAMHeader 
                  -> ConduitT (BAM, BAM) BED (ResourceT IO) ()
bamToUniqFragment tmpdir header = do
    tmp <- concatMapC toBed .| unlinesAsciiC .| sinkTempFile tmpdir "tmp.bed"
    tmpSort <- yield "" .| sinkTempFile tmpdir "tmp.bed"
    liftIO $ shelly $ escaping False $ run_ "sort" ["-T", T.pack tmpdir, "-k1,1"
            , T.pack tmp, ">", T.pack tmpSort]
    sourceFile tmpSort .| linesUnboundedAsciiC .| mapC f .|
        CL.groupBy ((==) `on` fst) .| mapC (\x -> score .~ Just (length x) $ snd $ head x)
  where
    f x = let (a, b) = B.break (=='\t') x
          in (a, fromLine $ B.tail b)
    toBed (a, b) = do
        chr1 <- refName header a
        chr2 <- refName header b
        let bed = toLine $ BED chr1 (min start1 start2) (max end1 end2) (Just $ queryName a) Nothing Nothing
        if chr1 == chr2
            then return $ key <> "\t" <> bed 
            else mzero
      where
        start1 = startLoc a
        end1 = endLoc a
        start2 = startLoc b
        end2 = endLoc b
        key = mkKey (a, b)
    mkKey x = let Pair{..} = makeKeyPair (const Nothing) x
                  show' y | y >= 0 = fromJust $ packDecimal y
                          | otherwise = "-" <> fromJust (packDecimal $ negate y)
              in B.concat
                    [ show' _ref_id1, "-"
                    , show' _ref_id2, "-"
                    , show' _loc1, "-"
                    , show' _loc2
                    , B.pack $ show _orientation
                    , if _leftmost then "1" else "0" ]
{-# INLINE bamToUniqFragment #-}