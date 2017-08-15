{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.Pipeline.Download where

import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import           Control.Monad
import qualified Data.ByteString.Char8 as B
import           Data.Coerce           (coerce)
import           Data.List.Split       (splitOn)
import           Data.Maybe
import           Data.Singletons       (SingI)
import qualified Data.Text             as T
import           Network.HTTP.Conduit
import           Shelly                hiding (FilePath)

sraToFastq :: SingI tags
           => FilePath
           -> File tags 'SRA
           -> IO (MaybePair '[Gzip] '[Gzip] 'Fastq)
sraToFastq outDir input = if input `hasTag` Pairend
    then Right <$> fastqDumpPair input
    else Left <$> fastqDump input
  where
    fastqDump fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ ".fastq.gz"
            f1 = location .~ f1_name $ fl
        case splitOn "+" (fl^.location) of
            [f] -> shelly $ run_ "fastq-dump" ["--origfmt", "-I", "--gzip",
                "-O" ,T.pack outDir, T.pack f]
            fs -> do
                forM_ fs $ \f -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                    "-O", T.pack outDir, T.pack f]
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ ".fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
        return (coerce f1 :: File '[Gzip] 'Fastq)
    fastqDumpPair fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ "_1.fastq.gz"
            f2_name = outDir ++ "/" ++ fl^.location ++ "_2.fastq.gz"
            f1 = location .~ f1_name $ fl
            f2 = location .~ f2_name $ fl

        case splitOn "+" (fl^.location) of
            [f] -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                "--split-files", "--gzip", "-O", T.pack outDir, T.pack f]
            fs -> do
                forM_ fs $ \f -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                    "--split-files", "-O", T.pack outDir, T.pack f]
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_1.fastq") fs
                    f2s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_2.fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
                    run_ "cat" $ f2s ++ ["|", "gzip", "-c", ">", T.pack f2_name]
                    run_ "rm" f2s
        return ( (coerce f1 :: File '[Gzip] 'Fastq)
               , (coerce f2 :: File '[Gzip] 'Fastq) )

               {-
downloadData :: FilePath
             -> SomeFile
             -> IO (Maybe SomeFile)
downloadData outDir (SomeFile fl) = case () of
    _ | getFileType fl == SRA -> if fl `hasTag` Pairend
            then Right <$> fastqDumpPair fl
            else Left <$> fastqDump fl
      | fl `hasTag` ENCODE -> do
  where
    downloadfile fl = do
      f1_name <- downloadENCODE (fl^.location) outDir
      let f1 = location .~ f1_name $ fl
      return $ if getFileType fl == Other
          then let ft = guessFormat f1_name
               in case toSing ft of
                    SomeSing ft' -> withSingI ft' $ SomeFile $ setFiletype ft' f1
          else SomeFile f1
    fastqDump fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ ".fastq.gz"
            f1 = location .~ f1_name $ fl
        case splitOn "+" (fl^.location) of
            [f] -> shelly $ run_ "fastq-dump" ["--origfmt", "-I", "--gzip",
                "-O" ,T.pack outDir, T.pack f]
            fs -> do
                forM_ fs $ \f -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                    "-O", T.pack outDir, T.pack f]
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ ".fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
        return $ SomeFile (coerce f1 :: File '[Gzip] 'Fastq)
    fastqDumpPair fl = do
        let f1_name = outDir ++ "/" ++ fl^.location ++ "_1.fastq.gz"
            f2_name = outDir ++ "/" ++ fl^.location ++ "_2.fastq.gz"
            f1 = location .~ f1_name $ fl
            f2 = location .~ f2_name $ fl

        case splitOn "+" (fl^.location) of
            [f] -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                "--split-files", "--gzip", "-O", T.pack outDir, T.pack f]
            fs -> do
                forM_ fs $ \f -> shelly $ run_ "fastq-dump" ["--origfmt", "-I",
                    "--split-files", "-O", T.pack outDir, T.pack f]
                let f1s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_1.fastq") fs
                    f2s = map (\x -> T.pack $ outDir ++ "/" ++ x ++ "_2.fastq") fs
                shelly $ escaping False $ do
                    run_ "cat" $ f1s ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                    run_ "rm" f1s
                    run_ "cat" $ f2s ++ ["|", "gzip", "-c", ">", T.pack f2_name]
                    run_ "rm" f2s
        return ( SomeFile (coerce f1 :: File '[Gzip] 'Fastq)
               , SomeFile (coerce f2 :: File '[Gzip] 'Fastq) )

    -}

-- | Download data from ENCODE portal to a given directory.
downloadENCODE :: String    -- ^ Accession number
               -> FilePath  -- ^ Output dir
               -> IO FilePath
downloadENCODE acc dir = do
     request <- parseRequest url
     manager <- newManager tlsManagerSettings
     runResourceT $ do
         response <- http request manager
         let filename = T.unpack $ snd $ T.breakOnEnd "filename=" $ T.pack $
                B.unpack $ fromJust $ lookup "Content-Disposition" $
                responseHeaders response
         responseBody response $$+- sinkFileBS (dir ++ "/" ++ filename)
         return $ dir ++ "/" ++ filename
  where
    url = "https://www.encodeproject.org/files/" ++ acc ++ "/@@download"
