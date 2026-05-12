-- | LEGACY / DEAD CODE — not compiled (this directory is not in the cabal
--   hs-source-dirs).
--
--   These are the original non-DirectSum implementations of `matrixPrep`'s
--   internals, plus a couple of DS helpers (`processStepDS`,
--   `updateInnerNodeVector`) that became unreachable when `createRotationsDS`
--   was refactored to use `foldBalanced` over a padded leaf array.
--
--   Kept here for reference in case anyone wants to revive the
--   Compose-and-control encoding (vs the DirectSum tree). To re-integrate:
--   - The non-DS path depends on `condMultQb`, `conditional`, `toBitString`
--     (still present in MatrixPreparation.hs as live helpers).
--   - `processStepDS` and `updateInnerNodeVector` reference the *original*
--     `createRotationsDS` foldl pattern; they would need rewiring against
--     today's `foldBalanced`-based version.
module Programs.MatrixPreparationLegacy where

import HQP
import Data.Complex
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.List
import Data.Maybe

import Numeric.LinearAlgebra (conj)
import HQP.QOp.MatrixSemantics (CMat)

import Programs.MatrixPreparation
  ( getRows, splitList, calculateComplexGate, calculateRealGate
  , toBitString, conditional, condMultQb
  )

----------------------------------------------------------------------
-- Non-DirectSum path (Compose + controlled gates instead of DirectSum tree)
----------------------------------------------------------------------

unitaryV :: CMat -> QOp
unitaryV mat =
    let rowsVec = getRows mat
        rowLen v = sqrt . V.sum $ V.map (\x -> x * conjugate x) v
        allLens = V.map rowLen rowsVec
        retQOp = buildRowQOp allLens
    in retQOp ⊗ I

unitaryU :: CMat -> QOp
unitaryU mat =
    let rowsVec  = getRows (conj mat)
        uBlocks  = V.map buildRowQOp rowsVec
        numQbits = ceiling (logBase 2 (fromIntegral (V.length uBlocks)) :: Double)
        indexedBlocks          = V.indexed uBlocks
        fixedUBlockProcessStep = uBlockProcessStep numQbits
        accQOp = V.foldl fixedUBlockProcessStep I indexedBlocks
        swap   = if numQbits == 1 then Permute [1,0]
                 else Permute ([numQbits .. 2*numQbits - 1] ++ [0 .. numQbits - 1])
    in accQOp <> swap

uBlockProcessStep :: Int -> QOp -> (Int, QOp) -> QOp
uBlockProcessStep numQbits accqOp (idx, qOp) =
    accqOp <> condMultQb (toBitString idx numQbits) conditional qOp

buildRowQOp :: Vector ComplexT -> QOp
buildRowQOp vs =
    let numQbits = ceiling (logBase 2 (fromIntegral (V.length vs)) :: Double)
        qOp      = Id numQbits
        pairList = splitList vs
    in createRotations numQbits 0 qOp pairList

createRotations :: Int -> Int -> QOp -> [V.Vector ComplexT] -> QOp
createRotations numQbits level inQOp vs
    | length vs == 1 = inQOp <> createQOp numQbits level 0 (vs !! 0)
    | otherwise =
        let (pairLenLst, accQOp) =
              foldl (\(accVec, qOp) (i, v) ->
                       processStep numQbits level (accVec, qOp) i v)
                    (V.empty, inQOp)
                    (zip [0..] vs)
        in createRotations numQbits (level + 1) accQOp (splitList pairLenLst)

processStep :: Int -> Int
            -> (Vector ComplexT, QOp) -> Int -> Vector ComplexT
            -> (Vector ComplexT, QOp)
processStep numQbits level (accVec, accQOp) pairIdx pairVector
    | V.length pairVector <= 2 = (newVec, newQOp)
    | otherwise = error "Her burde der kun være par eller singletons"
  where
    val0     = pairVector V.! 0
    val1     = fromMaybe 0 (pairVector V.!? 1)
    newValue = sqrt (val0 * conjugate val0 + val1 * conjugate val1)
    newVec   = V.snoc accVec newValue
    newQOp   = accQOp <> createQOp numQbits level pairIdx pairVector

createQOp :: Int -> Int -> Int -> V.Vector ComplexT -> QOp
createQOp numQbits level pairIdx pairVector
    | level >= numQbits = error "Ups ... level er ikke mindre end numQbits"
    | level < 0         = error "Ups ... level er negativt"
    | level == 0 =
        if abs r1 < 1e-9 && abs r2 < 1e-9 then I else
            let rotQOp = calculateComplexGate r1 phi1 r2 phi2
            in condMultQb (toBitString pairIdx (numQbits - level - 1)) conditional rotQOp
    | level > 0 =
        if abs r1 < 1e-9 && abs r2 < 1e-9 then I else
            let rotQOp   = calculateRealGate r1 r2
                condGate = condMultQb (toBitString (pairIdx `mod` 2 ^ (numQbits - level - 1)) (numQbits - level - 1)) conditional rotQOp
            in if numQbits - level == 1
                  then rotQOp ⊗ Id level
                  else condGate ⊗ Id level
    | otherwise = error "createQOp: unreachable"
  where
    val0       = pairVector V.! 0
    val1       = fromMaybe 0 (pairVector V.!? 1)
    (r1, phi1) = polar val0
    (r2, _)    = polar val1
    phi2       = snd (polar val1)

----------------------------------------------------------------------
-- DS helpers that became unreachable when createRotationsDS switched
-- from foldl to foldBalanced.
----------------------------------------------------------------------

-- | Original DS step from the foldl-based createRotationsDS. References the
-- (renamed) createQOp / left-leaning DirectSum, so does NOT compile against
-- the current MatrixPreparation.hs without rewiring.
--
-- processStepDS :: Int -> (Vector ComplexT, QOp) -> Vector ComplexT
--               -> (Vector ComplexT, QOp)
-- processStepDS level (accVec, accQOp) pairVector
--     | V.length pairVector <= 2 = (newVec, newQOp)
--     | otherwise = error "Her burde der kun være par eller singletons"
--   where
--     newVec = updateInnerNodeVector pairVector accVec
--     newQOp = DirectSum accQOp (createQOpDS level pairVector)

updateInnerNodeVector :: Vector ComplexT -> Vector ComplexT -> Vector ComplexT
updateInnerNodeVector pairVector accVec =
    let val0     = pairVector V.! 0
        val1     = fromMaybe 0 (pairVector V.!? 1)
        newValue = sqrt (val0 * conjugate val0 + val1 * conjugate val1)
    in V.snoc accVec newValue
