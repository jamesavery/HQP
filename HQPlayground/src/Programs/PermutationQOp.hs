module Programs.PermutationQOp where

import HQP

import qualified Data.Vector.Unboxed as V
import Data.Bit
import Data.Bits (xor)


------------------------------------------------------------------------------
-- Fra permutation til liste af transpositioner
-- Kommer her ... Tager en permuteret liste.
------------------------------------------------------------------------------


------------------------------------------------------------------------------
-- Fra transposition til QOp
------------------------------------------------------------------------------
-- Kommer her ...




-- Convert Nat to a vector of a specific length
toBitVecOfLength :: Nat -> Nat -> V.Vector Bit
toBitVecOfLength len n = V.generate len $ \i ->
    -- Extract the i-th bit of the number n
    if testBit n i then Bit True else Bit False

-- Main Logic
simulate :: Integer -> Integer -> IO ()
simulate nat1 nat2 = do
    -- Find the minimal bits needed to represent the larger number
    let bitLength = if nat1 == 0 && nat2 == 0 then 0 
                    else floor (logBase 2 (fromIntegral (max nat1 nat2))) + 1
    
    let v1 = toBitVecOfLength bitLength nat1
    let v2 = toBitVecOfLength bitLength nat2
    
    -- Pairwise modulo 2 (Length is preserved automatically)
    let result = zipBits xor v1 v2
    
    -- Step through the result (Length is still bitLength)
    V.forM_ result $ \b -> do
        -- Your simulation logic here
        print b
