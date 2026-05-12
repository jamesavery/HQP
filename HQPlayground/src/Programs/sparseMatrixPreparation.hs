module Programs.SparseMatrixPreparation where

import HQP
import HQP.QOp.MPSSemantics

import Programs.MatrixPreparation

import qualified Data.Vector as V
import GHC.TypeLits (Nat)




----------------------------------------------------------------------------------------
-- Unitary encoding U_A of 
----------------------------------------------------------------------------------------

sparseUAEncoding :: Nat -> QOp -> QOp -> QOp
sparseUAEncoding sQubits valOracle posOracle =
    let noQubits = op_qubits valOracle
        n = noQubits - (sQubits + 1)
        dOp = I ⊗ ( diffusionOp sQubits ) ⊗ (Id n)
    in 
       dOp ∘ (I ⊗ posOracle) ∘ valOracle ∘ dOp
    
diffusionOp :: Int -> QOp
diffusionOp n = foldr1 Tensor (replicate n H) 

----------------------------------------------------------------------------------------
-- positionOracle :: [QOp] -> QOp
-- This corresponds to the permutation operator O_c = ∑|j,c_jl> <l,j| 
-- (Theorem 4.1 Camps et al. arXiv:2203.10236v4)
-- Unitarity demands that the sparseColOracle c_jl is a permutation for each l
-- This implementation receives a vector of QOp permutations.
-- One permutation for each sparsity index l in [1,...,s]
-- The list is padded to 2^n = bit length of s
----------------------------------------------------------------------------------------
positionOracle :: [QOp] -> QOp
positionOracle perms = buildDirectSum perms

----------------------------------------------------------------------------------------
-- valueOracle :: GenericSparseMatrix -> QOp
----------------------------------------------------------------------------------------

valueOracle :: [QOp] -> QOp
valueOracle rotBlocks =
    let 
        -- Form the direct sum of matrices ... NAVNET PÅ METODEN ER SKIDT !!! Refakt
        uQop = runUntilFinal rotBlocks

        -- Construct a basis change permutation ... HMM har du styr på mængden af vektorer i basis her.....
        pQOp = Permute ([1,2^(n-1) + 1]++[2,2^(n-1) + 2]++[3,2^(n-1) + 3]++[4,2^(n-1) + 4])
    in 
        pQOp <> uQop <> (Transpose pQOp)

buildRotQOp :: ComplexT -> QOp
buildRotQOp val =
    let  
        (r1, phi1) = polar val
        r2 = sqrt(1 - (r1)^2)
        phi2 = 0
    in 
        calculateComplexGate r1 phi1 r2 phi2








----------------------------------------------------------------------------------------
-- Data structure for sparse matrix - Transposition type 
----------------------------------------------------------------------------------------

data GenericSparseMatrix = GenericSparseMatrix 
    { dataList   :: V.Vector ComplexT
    , rowIndices :: V.Vector Int
    , idxPtr     :: V.Vector Int 
    } deriving (Show)

-- Create an instance
myMatrix = GenericSparseMatrix [0.5, 0.2] [0, 1] [0, 1, 2] -- use vectors instead


-- For time being, all columns havee the same sparsity
getSparsity :: GenericSparseMatrix -> Nat
getSparsity s = (idxPtr mat V.! 1) - (idxPtr mat V.! 0)

getColQubits :: GenericSparseMatrix -> Nat
getColQubits mat = 
    let
        sparsity = getSparsity mat
    in
        ceiling $ (logBase 2 (fromIntegral (V.length (dataList mat)) / sparsity))


checkNoCollisions :: GenericSparseMatrix -> Bool
checkNoCollisions mat = all colUnique [0 .. V.length (idxPtr mat) - 2]
  where
    colUnique j = 
        let start = idxPtr mat V.! j
            end   = idxPtr mat V.! (j + 1)
            rows  = V.slice start (end - start) (rowIndices mat)
        -- Using V.uniq on a sorted slice is fastest, but nub is fine for small s
        in V.length (V.fromList . nub . V.toList $ rows) == V.length rows

checkLPermutation :: Int -> GenericSparseMatrix -> Int -> Bool
checkLPermutation n mat l = 
    let results = [getColOracle mat j l | j <- [0..n-1]]
        rowList = [r | Just r <- results]
    in length rowList == n && sort rowList == [0..n-1]

-- | Verify that every column has exactly s non-zero entries.
hasFixedSparsity :: Int -> GenericSparseMatrix -> Bool
hasFixedSparsity s mat = 
    -- We check every column index from 0 to N-1
    -- N is derived from the number of pointers (length idxPtr - 1)
    all (\j -> columnSize j == s) [0 .. V.length (idxPtr mat) - 2]
  where
    columnSize j = (idxPtr mat V.! (j + 1)) - (idxPtr mat V.! j)

-- | A "Smart Constructor" that validates all three conditions
validateMatrix :: Int -> Int -> GenericSparseMatrix -> Either String GenericSparseMatrix
validateMatrix n s mat

    | V.length (idxPtr mat) /= n + 1 = Left "Dimension mismatch: idxPtr size must be N + 1"
    | not (hasFixedSparsity s mat)   = Left $ "Sparsity mismatch: Not all columns have exactly " ++ show s ++ " elements"
    | not (hasNoRowCollisions mat)   = Left "Validation failed: Row collisions detected within a column"

    | not (all (checkLPermutation n mat) [0..s-1]) = Left "Validation failed: Fixed-l mapping is not a permutation"
    | otherwise = Right mat


----------------------------------------------------------------------------------------
-- sparseColOracle :: 
-- Returns the Row Index for the l-th non-zero element in column j
----------------------------------------------------------------------------------------

sparseColFunction :: GenericSparseMatrix -> Int -> Int -> Maybe Int
sparseColFunction mat l j = do
    startIdx <- idxPtr mat V.!? j
    endIdx   <- idxPtr mat V.!? (j + 1)
    let localOffset = startIdx + l
    if localOffset < endIdx
        then rowIndices mat V.!? localOffset
        else Nothing