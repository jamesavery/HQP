module Programs.SVT where
import HQP
import Programs.Qubitization
import Polynomial.Roots (roots)
import Polynomial.Basic (polyadd, polysub, polymult,polyderiv, polyeval)
import Data.Complex (Complex((:+)), realPart, imagPart,magnitude,polar)
import qualified Data.Vector as V
import Data.List (sortOn,sortBy, groupBy, partition)
import Data.Ord (comparing)
import Data.Sequence (Seq(..), (><))
import qualified Data.Sequence as Seq
import Debug.Trace

-- This vector ... 
svtVectorPhi :: [ComplexT] -> [ComplexT] -> Seq Double
svtVectorPhi p q =
    let
        -- 1. Generate the sequence
        -- (If qspVectorPhi returns a list, use Seq.fromList)
        qspVector = Seq.fromList (qspVectorPhi p q)

        -- 2. Modify to get the SVT vector
        (phi0, svtVector, phiLast) = case qspVector of
                    (p0 :<| (mid :|> pL))  -> (p0, mid, pL)
                    _                      -> error "Vector too short for SVT modification"

        d = (fromIntegral (Seq.length qspVector - 2) * (pi / 2) :: Double)
        modQspVector = fmap (subtract (pi / 2)) svtVector
    in
       Seq.fromList ([phi0 + phiLast + d ]) >< modQspVector

altPhaseMod :: Int -> Int -> Seq Double -> QOp -> QOp
altPhaseMod matrixQubits encQubits svtVec blockEnc =
    let
        -- Create the phase matrices
        phaseOperators = fmap (phaseOpBuilder matrixQubits encQubits) svtVec 
    in
        -- Create the phased alternating sequence
        phaseAltSeqBuilder (isEvenLength svtVec) phaseOperators blockEnc

phaseOpBuilder :: Int -> Int -> Double -> QOp
phaseOpBuilder matrixQubits encQubits phi = 
    let 
        phiPhaseOp = ((Id matrixQubits) ∘ Phase (toRational ( phi / pi))) ⊕ ((Id matrixQubits) ∘ Phase (toRational ( - phi / pi)))
    in
        phasePadder (matrixQubits + 1) (encQubits - 1) phi phiPhaseOp

phasePadder :: Int -> Int -> Double -> QOp -> QOp
phasePadder opQubits encQubits phi phiPhaseOp
    | encQubits == 0 = phiPhaseOp
    | otherwise =
        let 
            newphiPhaseOp = phiPhaseOp ⊕ (Id opQubits <> Phase (toRational ( - phi / pi)))
        in
            phasePadder (opQubits + 1) (encQubits - 1) phi newphiPhaseOp

phaseAltSeqBuilder :: Bool -> Seq QOp -> QOp -> QOp
phaseAltSeqBuilder isEven phaseOps blockEnc  
    | not isEven =  
        let
            (phaseOp, phaseOpsRest) = case phaseOps of
                (p0 :<| mid ) -> (p0, mid)
                _             -> error "Vector too short for SVT modification"

            part = phaseOp <> blockEnc
        in
            part <> phaseAltSeqBuilder True phaseOpsRest blockEnc

    | otherwise =
        let
            (phaseOp1, phaseOpsRest1) = case phaseOps of
                (p0 :<| rest ) -> (p0, rest)
                _              -> error "Vector too short for SVT modification"
            part1 = phaseOp1 <> (Adjoint blockEnc)

            (phaseOp2, phaseOpsRest2) = case phaseOpsRest1 of
                (p0 :<| rest ) -> (p0, rest)
                _              -> error "Vector too short for SVT modification"
            part2 = phaseOp2 <> blockEnc

            part = part1 <> part2
        
        in case phaseOpsRest2 of 
            Empty -> part 
            _  -> part <> phaseAltSeqBuilder True phaseOpsRest2 blockEnc
            

isEvenLength :: Seq a -> Bool
isEvenLength s = even (Seq.length s)