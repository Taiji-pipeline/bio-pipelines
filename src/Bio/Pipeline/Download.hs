{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Bio.Pipeline.Download
    ( downloadFiles
    , downloadENCODE
    , downloadENCODE'
    , sraToFastq
    , getUrl
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.File   (SFileType)
import           Bio.Data.Experiment.Parser (guessFormat)
import           Conduit
import           Control.Lens
import           Control.Monad
import qualified Data.ByteString.Char8      as B
import           Data.Coerce                (coerce)
import           Data.List                  (nub, isPrefixOf)
import           Data.List.Ordered          (nubSort)
import           Data.List.Split            (splitOn)
import Data.Conduit.Zlib (multiple, ungzip)
import Control.Concurrent.Async (forConcurrently_, forConcurrently, concurrently)
import           Data.Maybe
import Data.Char (toUpper)
import           Data.Singletons
import qualified Data.Text                  as T
import           Network.HTTP.Conduit
import           Shelly                     hiding (FilePath)
import Network.FTP.Client.Conduit (retr)
import Network.FTP.Client (withFTP, login)

import Bio.Pipeline.Utils (withTempDir)

downloadFiles :: FilePath
              -> FilePath
              -> Either SomeFile (SomeFile, SomeFile)
              -> IO (Either SomeFile (SomeFile, SomeFile))
downloadFiles outDir tempDir (Left (SomeFile fl))
    | getFileType fl == SRA = if fl `hasTag` PairedEnd
        then bimap SomeFile (bimap SomeFile SomeFile) <$>
            sraToFastq' outDir tempDir (coerce fl :: File '[PairedEnd] 'SRA)
        else bimap SomeFile (bimap SomeFile SomeFile) <$>
            sraToFastq' outDir tempDir (coerce fl :: File '[] 'SRA)
    | fl `hasTag` ENCODE = Left <$> downloadENCODE outDir (SomeFile fl)
    | map toUpper (take 4 $ fl^.location) == "URL:" = do
        let url = T.unpack $ T.strip $ T.pack $ drop 4 $ fl^.location
            output = outDir <> "/" <> T.unpack (snd $ T.breakOnEnd "/" $ T.pack $ fl^.location)
            newFile = location .~ output $ fl
        getUrl output url False
        return $ Left $ SomeFile newFile
    | otherwise = return $ Left $ SomeFile fl
downloadFiles outDir _ (Right (f1, f2)) = do
    f1' <- download f1
    f2' <- download f2
    return $ Right (f1', f2')
  where
    download (SomeFile fl)
        | fl `hasTag` ENCODE = downloadENCODE outDir (SomeFile fl)
        | map toUpper (take 4 $ fl^.location) == "URL:" = do
            let url = T.unpack $ T.strip $ T.pack $ drop 4 $ fl^.location
                output = outDir <> "/" <> T.unpack (snd $ T.breakOnEnd "/" $ T.pack $ fl^.location)
                newFile = location .~ output $ fl
            getUrl output url False
            return $ SomeFile newFile
        | otherwise = return $ SomeFile fl
 

downloadENCODE :: FilePath
               -> SomeFile
               -> IO SomeFile
downloadENCODE outDir (SomeFile (fl :: File filetag filetype)) = do
    let tags = fromSing (sing :: Sing filetag)
    f_name <- downloadENCODE' (fl^.location) outDir
    let newFile = location .~ f_name $ fl
        ft = guessFormat f_name
    return $ if getFileType fl == Other
        then case toSing ft of
                 SomeSing (ft' :: SFileType ft) -> withSingI ft' $ if gzipped f_name
                     then case toSing (nub $ Gzip : tags) of
                        SomeSing (tags' :: Sing filetag') -> withSingI tags' $
                            SomeFile (coerce newFile :: File filetag' ft)
                     else SomeFile (coerce newFile :: File filetag ft)
        else SomeFile newFile
  where
    gzipped = T.isSuffixOf ".gz" . T.pack
{-# INLINE downloadENCODE #-}

-- | Download data from ENCODE portal to a given directory.
downloadENCODE' :: String    -- ^ Accession number
                -> FilePath  -- ^ Output dir
                -> IO FilePath
downloadENCODE' acc dir = do
     request <- parseRequest url
     manager <- newManager tlsManagerSettings
     runResourceT $ do
         response <- http request manager
         let filename = T.unpack $ snd $ T.breakOnEnd "filename=" $ T.pack $
                B.unpack $ fromJust $ lookup "Content-Disposition" $
                responseHeaders response
         runConduit $ responseBody response .| sinkFileBS (dir ++ "/" ++ filename)
         return $ dir ++ "/" ++ filename
  where
    url = "https://www.encodeproject.org/files/" ++ acc ++ "/@@download"
{-# INLINE downloadENCODE' #-}

getUrl :: FilePath -> String -> Bool -> IO ()
getUrl output url inflate
    | "http" `isPrefixOf` url = httpDownload output url inflate
    | "ftp" `isPrefixOf` url = do
        let (_, url') = T.breakOnEnd "ftp://" (T.pack url)
            (base, filename) = T.breakOn "/" url'
        ftpDownload output (T.unpack base) (T.unpack filename) inflate
    | otherwise = error "Unknown protocol"

--------------------------------------------------------------------------------
-- Helper
--------------------------------------------------------------------------------

httpDownload :: FilePath
             -> String
             -> Bool
             -> IO ()
httpDownload output url inflate = do
    shelly $ mkdir_p $ fromText $ fst $ T.breakOnEnd "/" $ T.pack output
    request <- parseRequest url
    manager <- newManager tlsManagerSettings
    runResourceT $ do
        response <- http request manager
        runConduit $ responseBody response .| sink
  where
    sink | inflate = multiple ungzip .| sinkFileBS output
         | otherwise = sinkFileBS output
{-# INLINE httpDownload #-}

ftpDownload :: FilePath
            -> String
            -> String
            -> Bool
            -> IO ()
ftpDownload output base filename inflate = withFTP base 21 $ \h _ -> do
    _ <- login h "anonymous" ""
    shelly $ mkdir_p $ fromText $ fst $ T.breakOnEnd "/" $ T.pack output
    runResourceT $ runConduit $ retr h filename .| sink
  where
    sink | inflate = multiple ungzip .| sinkFileBS output
         | otherwise = sinkFileBS output
{-# INLINE ftpDownload #-}

sraToFastq :: (SingI tags, tags' ~ Insert' 'Gzip tags)
           => FilePath    -- ^ Output directory
           -> FilePath    -- ^ temp dir
           -> File tags 'SRA
           -> IO [File tags' 'Fastq]
sraToFastq outDir temp input = shelly (test_px "fasterq-dump") >>= \case
    False -> error "Please install the latest version of sra-tools: https://github.com/ncbi/sra-tools"
    True -> withTempDir (Just temp) $ \tmpDir -> do
        let tmpOutDir = tmpDir <> "/output/"
            sraAccessions = T.splitOn "+" $ T.pack $ input^.location
        forConcurrently_ sraAccessions $ \acc -> shelly $ run_ "fasterq-dump"
            [ "--split-files", "--include-technical"
            , "-O", T.pack tmpOutDir, acc, "-t", T.pack tmpDir ]
        suffixes <- nubSort . map (stripPrefixes sraAccessions . snd . T.breakOnEnd "/") <$>
            shelly (lsT tmpOutDir)
        forConcurrently suffixes $ \suffix -> do
            let files = map (\acc -> T.pack tmpOutDir <> acc <> suffix) sraAccessions
                output = outDir ++ "/" <> input^.location <> T.unpack suffix <> ".gz"
            shelly $ escaping False $ bashPipeFail bash_ "cat" $
                files ++ ["|", "gzip", "-c", ">", T.pack output]
            return $ coerce $ location .~ output $ input
  where
    stripPrefixes prefixes txt = fromMaybe (error $ "Unexpected file: " <> show txt) $
        foldl f Nothing prefixes
      where
        f Nothing prefix = T.stripPrefix prefix txt
        f suffix _ = suffix
{-# INLINE sraToFastq #-}

sraToFastq' :: SingI tags
            => FilePath    -- ^ Output directory
            -> FilePath    -- ^ temp dir
            -> File tags 'SRA
            -> IO (Either (File '[Gzip] 'Fastq)
                          (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
                  )
sraToFastq' outDir temp input = if input `hasTag` PairedEnd
    then Right <$> dumFqPair input
    else Left <$> dumpFq input
  where
    dumpFq fl = hasFasterqDump >>= \case
        True -> fasterqDump fl
        False -> fastqDump fl
    dumFqPair fl = hasFasterqDump >>= \case
        True -> fasterqDumpPair fl
        False -> fastqDumpPair fl
    hasFasterqDump = shelly $ test_px "fasterq-dump" >>= \case
        True -> return True
        False -> test_px "fastq-dump" >>= \case
            True -> return False
            False -> error "Please install sra-tools: https://github.com/ncbi/sra-tools"
    fasterqDump fl = withTempDir (Just temp) $ \tmpDir -> do
        let f1_name = outDir ++ "/" ++ fl^.location ++ ".fastq.gz"
            f1 = location .~ f1_name $ fl
        outputs <- forConcurrently (splitOn "+" $ fl^.location) $ \f -> do
            let output = T.pack $ outDir ++ "/" ++ f ++ ".fastq"
            shelly $ run_ "fasterq-dump"
                ["-O", T.pack outDir, T.pack f, "-t", T.pack tmpDir, "-f"]
            return output
        shelly $ escaping False $ do
            run_ "cat" $ outputs ++ ["|", "gzip", "-c", ">", T.pack f1_name]
            run_ "rm" outputs
        return (coerce f1 :: File '[Gzip] 'Fastq)
    fasterqDumpPair fl = withTempDir (Just temp) $ \tmpDir -> do
        let f1_name = outDir ++ "/" ++ fl^.location ++ "_1.fastq.gz"
            f2_name = outDir ++ "/" ++ fl^.location ++ "_2.fastq.gz"
        (outputs1, outputs2) <- fmap unzip $ forConcurrently (splitOn "+" $ fl^.location) $ \f -> do
            let f1 = T.pack $ outDir ++ "/" ++ f ++ "_1.fastq"
                f2 = T.pack $ outDir ++ "/" ++ f ++ "_2.fastq"
            shelly $ run_ "fasterq-dump"
                ["--split-files", "-O", T.pack outDir, T.pack f, "-t", T.pack tmpDir, "-f"]
            return (f1,f2)
        concurrently ( shelly $ escaping False $ do
                run_ "cat" $ outputs1 ++ ["|", "gzip", "-c", ">", T.pack f1_name]
                run_ "rm" outputs1 )
                ( shelly $ escaping False $ do
                    run_ "cat" $ outputs2 ++ ["|", "gzip", "-c", ">", T.pack f2_name]
                    run_ "rm" outputs2
                )
        return ( (coerce (location .~ f1_name $ fl) :: File '[Gzip] 'Fastq)
               , (coerce (location .~ f2_name $ fl) :: File '[Gzip] 'Fastq) )
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
{-# INLINE sraToFastq' #-}