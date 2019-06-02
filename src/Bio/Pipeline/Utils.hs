{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.Pipeline.Utils
    ( Directory
    , asDir
    , getPath
    , withTempDir
    , withTemp
    ) where

import           Control.Monad.IO.Class
import           Data.Aeson             (FromJSON (..), ToJSON (..))
import Data.Binary (Binary(..))
import           Data.String            (IsString (..))
import qualified Data.Text              as T
import           GHC.Generics           (Generic)
import           Shelly                 (fromText, mkdir_p, shelly)
import System.IO.Temp (withSystemTempDirectory, withTempDirectory, withSystemTempFile, withTempFile)
import Control.Monad.Catch (MonadMask)
import System.IO (hClose)

newtype Directory = Directory FilePath deriving (Show, Read, Ord, Eq, Generic)

instance IsString Directory where
    fromString = Directory

instance Semigroup Directory where
    (<>) (Directory x) (Directory y) = Directory (x ++ y)

instance ToJSON Directory where
    toEncoding (Directory fl) = toEncoding fl

instance FromJSON Directory where
    parseJSON = fmap Directory . parseJSON

instance Binary Directory

asDir :: FilePath -> Directory
asDir = Directory
{-# INLINE asDir #-}

-- | Return the path to the directory. This will create a new directory if it
-- doesn't exist.
getPath :: MonadIO m => Directory -> m FilePath
getPath (Directory dir) = liftIO $ shelly $
    mkdir_p (fromText $ T.pack dir) >> return dir
{-# INLINE getPath #-}

-- | Create a temporary directory in the given directory.
withTempDir :: (MonadMask m, MonadIO m)
            => Maybe FilePath   -- ^ Parent directory to create the directory in
            -> (FilePath -> m a) -> m a
withTempDir Nothing = withSystemTempDirectory "tmp_dir_"
withTempDir (Just dir) = withTempDirectory dir "tmp_dir_"
{-# INLINE withTempDir #-}

-- | Create a temporary file in the given directory.
withTemp :: (MonadMask m, MonadIO m)
         => Maybe FilePath   -- ^ Parent directory to create the directory in
         -> (FilePath -> m a) -> m a
withTemp Nothing f = withSystemTempFile "tmp_file_" $ \fl h ->
    liftIO (hClose h) >> f fl
withTemp (Just dir) f = withTempFile dir "tmp_file_" $ \fl h ->
    liftIO (hClose h) >> f fl
{-# INLINE withTemp #-}