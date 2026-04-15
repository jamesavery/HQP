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
--import qualified Data.IntMap as Seq

main :: IO ()
main = do
    --print (show (svtVectorPhiTest()))
    print (show (sqrtSVTtest()))
    print (show (evenPolySVTtest()))

svtVectorPhiTest :: () -> Seq Double
svtVectorPhiTest () =
    let 
        p = [0,0,0.5,0,0.5]
        (pC,qC) = qspPolys p
    in
        svtVectorPhi pC qC

sqrtSVTtest :: () -> [[ComplexT]]
sqrtSVTtest () =
    let 
        -- p(x) = x^2 giver ikke i sig selv de singulære værdier kvadreret ...men ...
        p = [0,0,1] -- x^2 ... som skal blive til pC = (1 + i)x^2 -i (even), qC = sqrt(2)ix (odd)
                    -- Hvilket er [0.0 :+ (-1.0),0.0 :+ 0.0,1.0 :+ 1.0]           ,[0.0 :+ 0.0,0.0 :+ sqrt(2)]
                    -- Jeg får:   [0.0 :+ (-1.0),0.0 :+ 0.0,1.0 :+ 1.0,0.0 :+ 0.0],[0.0 :+ 0.0,0.0 :+ 1.4142135623730951,0.0 :+ 0.0]
        (pC,qC) = qspPolys p
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

        -- ... tage gennemsnit af de tilstande man får ved anvendelse af 
        -- den alternerende fasemodulerende sekvens hørende til blockEnc med svtVec hhv (-svtVec) : 
        tempState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        tempState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1])))
        
        -- De singulære normeringskonstanter skal også kvadreres:
        finalState1 = ((frobeniusNormStrict matIn)^2 :+ 0) .* tempState1
        finalState2 = ((frobeniusNormStrict matIn)^2 :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        -- Det virker (sidst jeg prøvede :-) ... man får værdierne 4 og 1
        toLists matOut  

evenPolySVTtest :: () -> [[ComplexT]]
evenPolySVTtest () =
    let 
        -- Et lige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,0,0.5,0,0.5] -- 0.5x^2 + 0.5x^4 
        (pC,qC) = qspPolys p

        --p = [0.0 :+ (-1.0), 0.0 :+ 0.0.     , 1.0 :+ 1.0,0, 0, 0]
        --q = [0.0 :+ 0.0   , 0.0 :+ sqrt(2.0), 0.0 :+ 0.0]
        svtVec = svtVectorPhi pC qC

        -- Definer en matrix og sæt matrixQubits
        -- Dette er en ikke-diagonal matrix med singulære værdier 2 og 1
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

        -- ... tage gennemsnit af de tilstande man får ved anvendelse af 
        -- den alternerende fasemodulerende sekvens hørende til blockEnc med svtVec hhv (-svtVec) : 
        tempState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        tempState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1])))
        
        -- De singulære normeringskonstanter skal også kvadreres:
        normSV = 0.5*(frobeniusNormStrict matIn)^2 + 0.5*(frobeniusNormStrict matIn)^4

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        -- Virker det? Man burde få de singulære værdier 10 og 1 
        -- ... får 10.80000000 og 1.800000001 ... hvor komme de 0.8 fra?
        toLists matOut 






----------------------------- Helpers -----------------------------------------
-- Helper to round a single number to n decimal places
roundTo :: RealFloat a => Int -> a -> a
roundTo n x = fromIntegral (round (x * 10^n)) / (10^n)

-- Round both parts of a complex number
roundComplex :: RealFloat a => Int -> Complex a -> Complex a
roundComplex n (r :+ i) = (roundTo n r) :+ (roundTo n i)


