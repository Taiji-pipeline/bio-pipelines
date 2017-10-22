{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Bio.Pipeline.Download (downloadFiles) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.File   (SFileType)
import           Bio.Data.Experiment.Parser (guessFormat)
import           Conduit
import           Control.Lens
import           Control.Monad
import qualified Data.ByteString.Char8      as B
import           Data.Coerce                (coerce)
import           Data.List                  (nub)
import           Data.List.Split            (splitOn)
import           Data.Maybe
import           Data.Singletons
import qualified Data.Text                  as T
import           Network.HTTP.Conduit
import           Shelly                     hiding (FilePath)

sraToFastq :: SingI tags
           => FilePath    -- ^ Output directory
           -> File tags 'SRA
           -> IO (Either (File '[Gzip] 'Fastq)
                         (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
                 )
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
{-# INLINE sraToFastq #-}

downloadFiles :: FilePath
              -> Either SomeFile (SomeFile, SomeFile)
              -> IO (Either SomeFile (SomeFile, SomeFile))
downloadFiles outDir (Left (SomeFile fl))
    | getFileType fl == SRA = if fl `hasTag` Pairend
        then bimap SomeFile (bimap SomeFile SomeFile) <$>
            sraToFastq outDir (coerce fl :: File '[Pairend] 'SRA)
        else bimap SomeFile (bimap SomeFile SomeFile) <$>
            sraToFastq outDir (coerce fl :: File '[] 'SRA)
    | fl `hasTag` ENCODE = Left <$> downloadENCODE outDir (SomeFile fl)
    | otherwise = return $ Left $ SomeFile fl
downloadFiles outDir (Right (SomeFile f1, SomeFile f2))
    | f1 `hasTag` ENCODE && f2 `hasTag` ENCODE = do
        f1' <- downloadENCODE outDir (SomeFile f1)
        f2' <- downloadENCODE outDir (SomeFile f2)
        return $ Right (f1', f2')
    | otherwise = return $ Right (SomeFile f1, SomeFile f2)

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
    gzipped fl = ".gz" `T.isSuffixOf` T.pack fl
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
         responseBody response $$+- sinkFileBS (dir ++ "/" ++ filename)
         return $ dir ++ "/" ++ filename
  where
    url = "https://www.encodeproject.org/files/" ++ acc ++ "/@@download"
{-# INLINE downloadENCODE' #-}
