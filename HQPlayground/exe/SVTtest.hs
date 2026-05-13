{-# LANGUAGE OverloadedLists #-}
module Main where 

import HQP
import HQP.QOp.StatevectorSemantics

import Programs.Qubitization
import Programs.SVT
import Programs.MatrixPreparation

import Numeric.LinearAlgebra (fromLists, toLists, norm_Frob)

import Data.Complex
import Data.Sequence (Seq(..))

import Data.Poly

main :: IO ()
main = do
    print (show (svtVectorPhiTest()))
    print (show (sqrtSVTtest()))
    print (show (evenPolySVTtest()))
    print (show (oddPolySVTtest()))
    print (show (oddHermitianSVTtest()))
    print (show (evenHermitianSVTtest()))

svtVectorPhiTest :: () -> Seq Double
svtVectorPhiTest () =
    let 
        p = [0,0,0.5,0,0.5] -- p(x) = 1/2x^2+1/2x^4
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

        -- Her er en ikke-diagonal matrix med singulære værdier 2 og 1
        normMat = (1/(2 * sqrt(2)) :+ 0.0)
        matDataPre = [[(2.0 * sqrt(3.0) + 1.0)  :+ 0.0, (2.0 - sqrt(3.0)) :+ 0.0], [(2.0 * sqrt(3.0) - 1.0) :+ 0.0, (2.0 + sqrt(3.0)) :+ 0.0]]  
        matData = map (map (* normMat)) matDataPre

        matIn = fromLists matData
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
        finalState1 = ((norm_Frob matIn)^(2 :: Int) :+ 0) .* tempState1
        finalState2 = ((norm_Frob matIn)^(2 :: Int) :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
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
        normMat = (1/(2 * sqrt(2)) :+ 0.0)
        matDataPre = [[(2.0 * sqrt(3.0) + 1.0)  :+ 0.0, (2.0 - sqrt(3.0)) :+ 0.0], [(2.0 * sqrt(3.0) - 1.0) :+ 0.0, (2.0 + sqrt(3.0)) :+ 0.0]]  
        matData = map (map (* normMat)) matDataPre

        matIn = fromLists matData
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
        normSV = 0.5*(norm_Frob matIn)^(2 :: Int) + 0.5*(norm_Frob matIn)^(4 :: Int)

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        -- Virker det? Man burde få de singulære værdier 10 og 1 
        -- ... får 10.80000000 og 1.800000001 ... hvor komme de 0.8 fra?
        toLists matOut 

oddPolySVTtest :: () -> [[ComplexT]]
oddPolySVTtest () =
    let 
        -- Et ulige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,(-1.5),0,0.5,0,2] -- 2x^5 + 0.5x^3 - 1.5x 
        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Definer en matrix og sæt matrixQubits
        -- Dette er en ikke-diagonal matrix med singulære værdier 0.8 og 0.6
        matDataPre =   [[((4*sqrt(6) - 3*sqrt(2))/20) :+ 0.0, ((-4*sqrt(6)-3*sqrt(2))/20) :+ 0.0], 
                        [((4*sqrt(2)+3*sqrt(6))/20) :+ 0.0, ((-4*sqrt(2)+3*sqrt(6))/20) :+ 0.0]]  
        matData = matDataPre

        matIn = fromLists matData
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
        
        -- De singulære normeringskonstanter beregnes:
        normSV = 0.5*(norm_Frob matIn)^(2 :: Int) + 0.5*(norm_Frob matIn)^(4 :: Int)

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        -- Virker det? Man burde få de singulære værdier 0.5248 og 0.2448. Det får man 
        toLists matOut 

oddHermitianSVTtest :: () -> [[ComplexT]]
oddHermitianSVTtest () =
    let 
        -- Et ulige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,(-1.5),0,0.5,0,2] -- 2x^5 + 0.5x^3 - 1.5x 
        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Definer en matrix og sæt matrixQubits
        matDataPre =   [[0.0 :+ 0.0, sqrt(2)/4 :+ (-sqrt(2)/4)], 
                        [sqrt(2)/4 :+ sqrt(2)/4, (-1/sqrt(2)) :+ 0.0]]  
        matData = matDataPre

        matIn = fromLists matData
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
        
        -- De singulære normeringskonstanter beregnes:
        normSV = 0.5*(norm_Frob matIn)^(2 :: Int) + 0.5*(norm_Frob matIn)^(4 :: Int)

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        -- Virker det? Man burde få de singulære værdier 0.5248 og 0.2448. Det får man 
        toLists matOut 

evenHermitianSVTtest :: () -> [[ComplexT]]
evenHermitianSVTtest () =
    let 
        -- Et lige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,0,0.5,0,0.5] :: VPoly Rational -- 0.5x^2 + 0.5x^4 

        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Definer en Hermitisk matrix og sæt matrixQubits
        matDataPre =   [[0.0 :+ 0.0, sqrt(2)/4 :+ (-sqrt(2)/4)], 
                        [sqrt(2)/4 :+ sqrt(2)/4, (-1/sqrt(2)) :+ 0.0]]   
        matData = matDataPre

        matIn = fromLists matData
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
        
        -- De singulære normeringskonstanter beregnes:
        normSV = 0.5*(norm_Frob matIn)^(2 :: Int) + 0.5*(norm_Frob matIn)^(4 :: Int)

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData 
    in 
        -- Virker det? Man burde få de singulære værdier 0.5248 og 0.2448. Det får man 
        toLists matOut 




----------------------------- Helpers -----------------------------------------
-- Helper to round a single number to n decimal places
roundTo :: (RealFloat a) => Int -> a -> a
roundTo n x = fromIntegral (round (x * 10^n) :: Integer) / (10^n)

roundComplex :: (RealFloat a) => Int -> Complex a -> Complex a
roundComplex n (r :+ i) = (roundTo n r) :+ (roundTo n i)


