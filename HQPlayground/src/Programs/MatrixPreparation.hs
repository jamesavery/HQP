module Programs.MatrixPreparation where
import HQP
import Data.Complex
import Data.Ratio
import System.Random(mkStdGen, randoms)
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.List
import Data.Maybe
import Debug.Trace

import Numeric.LinearAlgebra (conj, toRows, toList)
import HQP.QOp.MatrixSemantics (CMat)

-- | The row-prep internals below consume Data.Vector, so convert from hmatrix's
--   storable rows once here. Remove when the row-prep is refactored too.
getRows :: CMat -> V.Vector (V.Vector ComplexT)
getRows = V.fromList . map (V.fromList . toList) . toRows

----------------------Matrixpreparation--------------------------

unitaryV :: CMat -> QOp
unitaryV mat = 
    let rowsVec = getRows mat
        rowLen v = sqrt . V.sum $ V.map (\x -> x * conjugate x) v
        allLens = V.map rowLen rowsVec
        retQOp = buildRowQOp allLens
    in 
        retQOp ⊗ I

-- Direct sum version
unitaryVDS :: CMat -> QOp
unitaryVDS mat = 
    let rowsVec = getRows mat
        rowLen v = sqrt . V.sum $ V.map (\x -> x * conjugate x) v
        allLens = V.map rowLen rowsVec
        retQOp = buildRowQOpDS allLens
        numQbits = ceiling (logBase 2 (fromIntegral (V.length rowsVec)))
    in 
        retQOp ⊗ (Id numQbits)

unitaryU:: CMat -> QOp
unitaryU mat =  
    let
        rowsVec = getRows (conj mat)
        uBlocks = (V.map buildRowQOp rowsVec)
        numQbits = ceiling (logBase 2 (fromIntegral (V.length uBlocks)))
        indexedBlocks = V.indexed uBlocks
        fixedUBlockProcessStep = uBlockProcessStep numQbits
        accQOp = V.foldl fixedUBlockProcessStep I indexedBlocks 
        swap = if numQbits == 1 then Permute [1,0] else Permute ([numQbits .. 2*(numQbits) - 1 ] ++ [0 .. (numQbits - 1)])
    in
        accQOp <> swap


uBlockProcessStep :: Int -> QOp -> (Int,QOp) -> QOp
uBlockProcessStep numQbits accqOp (idx, qOp) = 
        accqOp <> condMultQb (toBitString idx numQbits) conditional qOp


-- Direct sum version
unitaryUDS:: CMat -> QOp
unitaryUDS mat =  
    let
        -- Get the row encodings
        rowsVec = getRows (conj mat)
        uBlocks = (V.map buildRowQOpDS rowsVec)
        -- Pad with I_n to 2^n length
        numQbits = ceiling (logBase 2 (fromIntegral (V.length uBlocks)))
        uBlocksPadded = padToPowerOf2 numQbits (Id numQbits) uBlocks
        -- Build the direct sums
        dsQOp = runUntilFinal uBlocksPadded
        --dsQOp = V.foldl1 processPair uBlocksPadded

        -- Building the SWAP
        swap = if numQbits == 1 then Permute [1,0] else Permute ([numQbits .. 2*(numQbits) - 1 ] ++ [0 .. (numQbits - 1)])
    in
        dsQOp <> swap

padToPowerOf2 :: Int -> QOp -> V.Vector QOp -> V.Vector QOp
padToPowerOf2 numQbits paddingObj vec
    | len == (2^numQbits) = vec -- Already 2^n
    | otherwise   = vec V.++ V.replicate ((2^numQbits) - len) paddingObj
  where
    len  = V.length vec

-- 1. Operation for intermediate pairs
processPair :: QOp -> QOp -> QOp
processPair x y = DirectSum x y

reduceStep :: V.Vector QOp -> V.Vector QOp
reduceStep vec = V.generate halfSize $ \i ->
    let op1 = vec V.! (2 * i)
        op2 = vec V.! (2 * i + 1)
    in processPair op1 op2
  where
    halfSize = V.length vec `div` 2

-- 3. The recursive driver
runUntilFinal :: V.Vector QOp -> QOp
runUntilFinal vec
    | V.length vec == 2 = processPair (vec V.! 0) (vec V.! 1) 
    | otherwise         = runUntilFinal (reduceStep vec)


matrixPrep :: CMat -> QOp
matrixPrep mat = 
    let
        -- Find U matrix
        --uQOp = unitaryU mat
        uQOp = unitaryUDS mat -- Direct Sum version

        -- Find V matrix
        --vQOp = unitaryV mat
        vQOp = unitaryVDS mat
    in
        (Adjoint uQOp) ∘ vQOp

------------------------ Row preparation -------------------------

buildRowQOp :: Vector ComplexT -> QOp
buildRowQOp vs = do
    let
        numQbits = ceiling (logBase 2 (fromIntegral (V.length vs)))
        qOp = Id (numQbits)   
        pairList = splitList vs
    
    createRotations numQbits 0 qOp pairList

-- Direct sum version
buildRowQOpDS :: Vector ComplexT -> QOp
buildRowQOpDS vs = do
    let pairList = splitList vs
    createRotationsDS 0 pairList

createRotations :: Int -> Int -> QOp -> [V.Vector ComplexT] -> QOp
createRotations numQbits level inQOp vs

    | length vs == 1 = inQOp <> createQOp numQbits level 0 (vs !! 0)
    | otherwise =
        let (pairLenLst, accQOp) = foldl (\(accVec, qOp) (i, v) -> 
                                    processStep numQbits level (accVec, qOp) i v) 
                                   (V.empty, inQOp) 
                                   (zip [0..] vs)
        in createRotations numQbits (level + 1) accQOp (splitList pairLenLst)

-- Direct sum version
createRotationsDS :: Int -> [V.Vector ComplexT] -> QOp
createRotationsDS level vs
    | length vs == 1 = (createQOpDS level (head vs) ⊗ (Id level))
    | otherwise =
        let (pairLenLst, accQOp) = case uncons vs of
                Just (first, rest) -> 
                    -- Initialize with the first element
                    let iniVec = updateInnerNodeVector first V.empty
                        iniQOp = createQOpDS level first
                    in foldl (\(accVec, qOp) v -> 
                                processStepDS level (accVec, qOp) v) 
                             (iniVec, iniQOp) 
                             rest

                Nothing -> error "createRotationsDS: list was unexpectedly empty."
        
        -- Recursive call and composition
        in (accQOp ⊗ (Id level)) <> createRotationsDS (level + 1) (splitList pairLenLst)  

-- Direct Sum Version
processStepDS :: Int -> (Vector ComplexT, QOp) -> Vector ComplexT -> (Vector ComplexT, QOp)
processStepDS level (accVec, accQOp) pairVector
    | V.length pairVector <= 2 = 
        (newVec, newQOp) 
    | otherwise = error "Her burde der kun være par eller singletons"
  where
    newVec   = updateInnerNodeVector pairVector accVec 
    newQOp   = DirectSum accQOp (createQOpDS level pairVector)

updateInnerNodeVector:: Vector ComplexT -> Vector ComplexT -> Vector ComplexT
updateInnerNodeVector pairVector accVec =
    let
        val0     = pairVector V.! 0
        val1     = fromMaybe 0 (pairVector V.!? 1)
        newValue = sqrt (val0 * (conjugate val0) + val1 * (conjugate val1))
    in V.snoc accVec newValue

processStep :: Int -> Int -> (Vector ComplexT, QOp) -> Int -> Vector ComplexT -> (Vector ComplexT, QOp)
processStep numQbits level (accVec, accQOp) pairIdx pairVector
    | V.length pairVector <= 2 = 
        (newVec, newQOp) 
    | otherwise = error "Her burde der kun være par eller singletons"
  where
    val0     = pairVector V.! 0
    val1     = fromMaybe 0 (pairVector V.!? 1)
    newValue = sqrt (val0 * (conjugate val0) + val1 * (conjugate val1))
    newVec   = V.snoc accVec newValue
    newQOp   = accQOp <> createQOp numQbits level pairIdx pairVector



createQOp :: Int -> Int -> Int -> V.Vector ComplexT -> QOp
createQOp numQbits level pairIdx pairVector
    -- 1. Error checks first

    | level >= numQbits = error "Ups ... level er ikke mindre end numQbits"
    | level < 0         = error "Ups ... level er negativt"
    
    -- 2. Base case

    | level == 0 = 
        if abs r1 < 1e-9 && abs r2 < 1e-9 then I else
            let rotQOp = calculateComplexGate r1 phi1 r2 phi2
            in condMultQb (toBitString pairIdx (numQbits - level -1)) conditional rotQOp

    -- 3. Recursive/Higher levels
    | level > 0 = 
        if abs r1 < 1e-9 && abs r2 < 1e-9 then I else
            let rotQOp = calculateRealGate r1 r2
                condGate = condMultQb (toBitString (pairIdx `mod` 2^(numQbits - level - 1)) (numQbits-level - 1)) conditional rotQOp
            in if (numQbits - level == 1) 
            then rotQOp ⊗ (Id level) 
            else condGate ⊗ (Id level)
    where
        val0       = pairVector V.! 0
        val1       = fromMaybe 0 (pairVector V.!? 1) 
        (r1, phi1) = polar val0
        (r2, phi2) = polar val1

-- Direct Sum Version
createQOpDS :: Int -> V.Vector ComplexT -> QOp
createQOpDS level pairVector
    | level < 0         = error "Ups ... level er negativt"
    | level == 0 = 
        if abs r1 < 1e-9 && abs r2 < 1e-9 then I else calculateComplexGate r1 phi1 r2 phi2
    | level > 0 = 
        if abs r1 < 1e-9 && abs r2 < 1e-9 then I else calculateRealGate r1 r2
    where
        val0       = pairVector V.! 0
        val1       = fromMaybe 0 (pairVector V.!? 1) 
        (r1, phi1) = polar val0
        (r2, phi2) = polar val1


calculateComplexGate :: Double -> Double -> Double -> Double -> QOp
calculateComplexGate r1 phi1 r2 phi2 
    | abs r1 < 1e-9 && abs r2 < 1e-9 = I
    | otherwise =
        let phi    = toRational $   (phi2 - phi1)/pi
            lambda = toRational $  -(phi1 + phi2)/pi
            theta  = toRational $ if r1 == 0 then 1 else (2/pi) * acos (r1 / sqrt (r1**2 + r2**2)) 
        in 
            (R Z (-phi)) ∘ (R Y (-theta)) ∘ (R Z (-lambda))

calculateRealGate :: Double -> Double -> QOp
calculateRealGate r1 r2
    | abs r1 < 1e-9 && abs r2 < 1e-9 = I
    | otherwise =
        let theta = toRational $ if abs r1 < 1e-9 then 1 else (2/pi) * acos (r1 / sqrt (r1**2 + r2**2)) 
        in 
            R Y (-theta)

toBitString :: Int -> Int -> [Int]
toBitString n bitStrLen

    | n < 0     = replicate bitStrLen 0  -- Or handle error as needed
    | otherwise = 
        let 
            -- Generate bits from least to most significant
            bits = unfoldr step n
            step 0 = Nothing
            step x = Just (fromIntegral (x `mod` 2), x `div` 2)
            
            -- Reverse to get proper order and pad with leading zeros
            rawBits = reverse bits
            padding = replicate (bitStrLen - length rawBits) 0
            
            -- Combine and ensure final length is exactly bitStrLen
        in take bitStrLen (padding ++ rawBits)

conditional :: Int -> QOp -> QOp
conditional bit inQOp 
    | bit == 0 = (X ⊗ I) <> C inQOp <> (X ⊗ I)
    | bit == 1 = C inQOp
    | otherwise = error "Ups ... bit skal være 0 eller 1" 

condMultQb :: [Int] -> (Int -> QOp -> QOp) -> QOp -> QOp
condMultQb bitPattern fct initialOp =
    foldr (\bit op -> fct bit op) initialOp bitPattern

----------------- Helper funcs ----------
-- Burde nok slås sammen med koden ved runUntilFinal
-- Det er en træopbyggende rekursiv metode der bør abstraheres

splitList :: V.Vector ComplexT -> [V.Vector ComplexT]
splitList vs = go [vs]
  where
    go curr

        | all (\v -> V.length v <= 2) curr = curr
        | otherwise = go (concatMap splitStep curr)

-- 2. Logic to split a single Vector based on the 2^n rule
splitStep :: V.Vector ComplexT -> [V.Vector ComplexT]
splitStep v

    | len <= 2 = [v]
    | otherwise =
        let n = floor (logBase 2 (fromIntegral len) :: Double)
            -- If len is exactly a power of 2 (e.g., 4), split it in half (2, 2)
            -- Otherwise, split at the largest power of 2 (e.g., 6 becomes 4, 2)
            splitPoint = if 2 ^ n == len then 2 ^ (n - 1) else 2 ^ n
            left  = V.slice 0 splitPoint v
            right = V.slice splitPoint (len - splitPoint) v
        in [left, right]
  where len = V.length v







