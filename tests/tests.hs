import Test.Tasty
import Test

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ tests ]
