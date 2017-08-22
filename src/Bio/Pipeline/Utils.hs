{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.Pipeline.Utils
    ( Directory
    , asDir
    , getPath
    ) where

import           Control.Monad.IO.Class
import           Data.Aeson (ToJSON(..), FromJSON(..))
import           Data.Serialize         (Serialize (..))
import qualified Data.Text              as T
import           GHC.Generics           (Generic)
import           Shelly                 (fromText, mkdir_p, shelly)

newtype Directory = Directory FilePath deriving (Show, Read, Ord, Eq, Generic)

instance ToJSON Directory where
    toEncoding (Directory fl) = toEncoding fl

instance FromJSON Directory where
    parseJSON = fmap Directory . parseJSON

instance Serialize Directory

asDir :: FilePath -> Directory
asDir = Directory

-- | Return the path to the directory. This will create a new directory if it
-- doesn't exist.
getPath :: MonadIO m => Directory -> m FilePath
getPath (Directory dir) = liftIO $ shelly $
    mkdir_p (fromText $ T.pack dir) >> return dir
