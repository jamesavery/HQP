module Main where 

import HQP
import HQP.QOp
import HQP.QOp.StatevectorSemantics
import HQP.PrettyPrint

import Programs.Qubitization
import Programs.SVT
import Programs.MatrixPreparation

import Data.Complex (Complex((:+)),magnitude)
import Polynomial.Roots (roots) -- From dsp package
import Polynomial.Basic (polyadd, polysub, polymult)
import Data.Sequence (Seq(..), (><))
import qualified Data.Sequence as Seq

import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.Complex
import Data.List
import Data.Maybe

import Debug.Trace

main :: IO ()
main = do
    --print (show (svtVectorPhiTest()))
    print (show (altPhaseModTest()))


svtVectorPhiTest :: () -> Seq Double
svtVectorPhiTest () =
    let 
        p = [0.0 :+ (-1.0), 0, 1.0 :+ 1.0,0, 0, 0]
        q = [0, 0.0 :+ sqrt(2.0),0]
    in
        svtVectorPhi p q


altPhaseModTest :: () -> [[ComplexT]]
altPhaseModTest () =
    let 
        -- p(x) = x^2 giver ikke i sig selv de singulære værdier kvadreret ...men ...
        p = [0,0,1] -- x^2 ... som skal blive til pC = (1 + i)x^2 -i (even), qC = sqrt(2)ix (odd)
                    -- Hvilket er [0.0 :+ (-1.0),0.0 :+ 0.0,1.0 :+ 1.0]           ,[0.0 :+ 0.0,0.0 :+ sqrt(2)]
                    -- Jeg får:   [0.0 :+ (-1.0),0.0 :+ 0.0,1.0 :+ 1.0,0.0 :+ 0.0],[0.0 :+ 0.0,0.0 :+ 1.4142135623730951,0.0 :+ 0.0]
        (pC,qC) = qspPolys p

        --p = [0.0 :+ (-1.0), 0.0 :+ 0.0.     , 1.0 :+ 1.0,0, 0, 0]
        --q = [0.0 :+ 0.0   , 0.0 :+ sqrt(2.0), 0.0 :+ 0.0]
        svtVec = svtVectorPhi pC qC

        -- Definer en matrix og sæt matrixQubits
        -- Her er en diagonalmatrix med de singulære værdier 0.6 og 0.3
        -- (Testen virker her ... men det er måske ikke så imponerende :-)
        --matData =  [[ 0.6 :+ 0.0, 0.0 :+ 0.0 ], 
        --            [ 0.0 :+ 0.0, 0.3 :+ 0.0 ]]

        -- Her er en ikke-diagonal matrix med singulære værdier 2 og 1
        norm = (1/(2 * sqrt(2)) :+ 0.0)
        matDataPre = [[(2.0 * sqrt(3.0) + 1.0)  :+ 0.0, (2.0 - sqrt(3.0)) :+ 0.0], [(2.0 * sqrt(3.0) - 1.0) :+ 0.0, (2.0 + sqrt(3.0)) :+ 0.0]]  
        matData = map (map (* norm)) matDataPre

        matIn = fromRows matData
        matrixQubits = 1

        -- Lav block encoding og sæt encQubits
        blockEnc = matrixPrep matIn
        encQubits = 1

        -- SVT udføres
        svtQOp = altPhaseMod matrixQubits encQubits svtVec blockEnc
        svtQOt = evalOp $ svtQOp

        -- SVT udføres med (-1 * svtVec) ... med det formål at ...
        svtVecMinus = fmap (* (-1)) svtVec
        svtConjQOp = altPhaseMod matrixQubits encQubits svtVecMinus blockEnc
        svtConjQOt = evalOp $ svtConjQOp

        -- ... tage gennemsnit af de to tilstande hørende til 
        -- svtVec og (-svtVec) : 
        tempState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        tempState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1])))
        
        -- De singulære normeringskonstanten skal også kvadreres:
        finalState1 = ((frobeniusNormStrict matIn)^2 :+ 0) .* tempState1
        finalState2 = ((frobeniusNormStrict matIn)^2 :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2

        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix har de singulære værdier i diagonalen
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        -- Det virker ... jeg får værdierne 4 og 1
        toLists matOut  

-- Helper to round a single number to n decimal places
roundTo :: RealFloat a => Int -> a -> a
roundTo n x = fromIntegral (round (x * 10^n)) / (10^n)

-- Round both parts of a complex number
roundComplex :: RealFloat a => Int -> Complex a -> Complex a
roundComplex n (r :+ i) = (roundTo n r) :+ (roundTo n i)


