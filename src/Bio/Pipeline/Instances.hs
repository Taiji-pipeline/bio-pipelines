{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TemplateHaskell    #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

module Bio.Pipeline.Instances where

import           Bio.Data.Bed.Types
import           Bio.Data.Bed.Utils
import           Bio.Motif             (Motif (..), PWM (..), Bkgd(..), CDF(..))
import           Control.Monad
import           Data.Aeson
import           Data.Aeson.TH
import           Data.Binary           (Binary)
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk, original)
import qualified Data.Matrix.Unboxed   as MU
import           Data.Serialize        (Serialize (..))
import qualified Data.Text             as T
import           Data.Vector.Binary    ()
import           Data.Vector.Serialize ()
import           GHC.Generics          (Generic)

instance FromJSON B.ByteString where
    parseJSON (String x) = return $ B.pack $ T.unpack x
    parseJSON _          = mzero

instance ToJSON B.ByteString where
    toJSON x = String $ T.pack $ B.unpack x

instance FromJSON (CI B.ByteString) where
    parseJSON = fmap mk . parseJSON

instance ToJSONKey (CI B.ByteString)

instance FromJSONKey (CI B.ByteString)

instance ToJSON (CI B.ByteString) where
    toJSON = toJSON . original

instance FromJSON (MU.Matrix Double) where
    parseJSON = genericParseJSON defaultOptions

instance ToJSON (MU.Matrix Double) where
    toEncoding = genericToEncoding defaultOptions

$(deriveJSON defaultOptions ''BED)
$(deriveJSON defaultOptions ''BED3)
$(deriveJSON defaultOptions ''NarrowPeak)
$(deriveJSON defaultOptions ''PWM)
$(deriveJSON defaultOptions ''Motif)

deriving instance Generic BED
instance Serialize BED
instance Binary BED

deriving instance Generic (BEDExt a b)
instance (Serialize a, Serialize b) => Serialize (BEDExt a b)
instance (Binary a, Binary b) => Binary (BEDExt a b)

instance Serialize (MU.Matrix Double)
instance Binary (MU.Matrix Double)

deriving instance Generic PWM
instance Serialize PWM
instance Binary PWM

deriving instance Generic Motif
instance Serialize Motif
instance Binary Motif

deriving instance Generic Bkgd
instance Serialize Bkgd
instance Binary Bkgd

deriving instance Generic CDF
instance Serialize CDF
instance Binary CDF

deriving instance Generic CutoffMotif
instance Serialize CutoffMotif
instance Binary CutoffMotif

instance Serialize (CI B.ByteString) where
    put = put . original
    get = fmap mk get
