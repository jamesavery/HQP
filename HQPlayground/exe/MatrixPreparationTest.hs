module Main where 

import HQP.QOp
import HQP.QOp.MPSSemantics
import HQP.PrettyPrint

import Programs.MatrixPreparation
import Programs.MatrixArithmetic

import Control.Monad (forM_)

import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.Complex
import Data.List
import Data.Maybe

import Debug.Trace


printbuildRowQOpTester :: IO ()
printbuildRowQOpTester = do
    let results = buildRowQOpTester ()
    forM_ results $ \res -> do
        printS res

main :: IO ()
main = do
    --printS $ ( (matrixPrepTester2 ())) 
    --print (show (matrixPrepTester3 ())) 
    --print (show (matrixPrepTester4 ())) 
    --printS $ ( (matrixPrepTester4 ()))
    --print (show (matrixPrepTester5 ())) 
    --printS $ ( (matrixPrepTester5 ()))
    printS $ ( (productTester ()))
    --print (show (productTester ())) --DENNE GIVER MÆRKELIGT OUTPUT. VIS FREM
    --print (show (sumTester ()))
    --printbuildRowQOpTester
    --printS $ (directSumTest ()) -- DENNE GIVER VIST OGSÅ MÆRKELIGT OUTPUT
 
 {--
    let 
        opInner001 = ((X ⊗ I) <> (C I) <> (X ⊗ I))
        op001 = (C opInner001)
        op001Alt = (C (Id 2))  

    printS $ ((apply (evalOp opInner001) (ket [0,0])))
    printS $ ((apply (evalOp opInner001) (ket [0,1])))
    printS $ ((apply (evalOp opInner001) (ket [1,0])))
    printS $ ((apply (evalOp opInner001) (ket [1,1])))
    print "---------Her kommer et underligt resultat:-----------"
    printS $ ((apply (evalOp op001) (ket [0,0,0])))
    printS $ ((apply (evalOp op001) (ket [0,0,1])))
    printS $ ((apply (evalOp op001) (ket [0,1,0])))
    printS $ ((apply (evalOp op001) (ket [0,1,1])))
    printS $ ((apply (evalOp op001) (ket [1,0,0])))
    printS $ ((apply (evalOp op001) (ket [1,0,1])))
    printS $ ((apply (evalOp op001) (ket [1,1,0])))
    printS $ ((apply (evalOp op001) (ket [1,1,1])))
    print "---------Her er hvordan det burde se ud:-----------"
    printS $ ((apply (evalOp op001Alt) (ket [0,0,0])))
    printS $ ((apply (evalOp op001Alt) (ket [0,0,1])))
    printS $ ((apply (evalOp op001Alt) (ket [0,1,0])))
    printS $ ((apply (evalOp op001Alt) (ket [0,1,1])))
    printS $ ((apply (evalOp op001Alt) (ket [1,0,0])))
    printS $ ((apply (evalOp op001Alt) (ket [1,0,1])))
    printS $ ((apply (evalOp op001Alt) (ket [1,1,0])))
    printS $ ((apply (evalOp op001Alt) (ket [1,1,1])))
--}

directSumTest :: () -> StateT
directSumTest () =
    let vec = V.fromList [X,X,X,X]
        dsQOp =  V.foldl1 processPairTest vec
        dsQOt = trace(show dsQOp) $ evalOp $ dsQOp 
    in
        apply dsQOt (ket [1,0,1,0])

processPairTest :: QOp -> QOp -> QOp
processPairTest x y = DirectSum x y


sumTester :: () -> [[ComplexT]]
sumTester () =
    let
        -- Matrix data
        matData1 =   [[1,2,3], [4,5,6], [7,8,9]]
        matData2 =   [[10,11,12], [13,14,15],[16,17,18]]    -- sum = [[11,13,15], [17,19,21],[23,25,27]]
        matIn1 = fromRows matData1
        matIn2 = fromRows matData2
        
        -- Norms
        norm1 = (frobeniusNormStrict matIn1 :+ 0)
        norm2 = (frobeniusNormStrict matIn2 :+ 0)

        -- Norm QOp
        normQOp = buildRowQOpDS (V.fromList [sqrt(norm1), sqrt(norm2)])  

        matData1norm = map (map (/ norm1)) [[1,2,3], [4,5,6], [7,8,9]]
        matData2norm = map (map (/ norm2)) [[10,11,12], [13,14,15],[16,17,18]]
        
        -- Normerede data matricer
        matIn1norm = fromRows matData1norm
        matIn2norm = fromRows matData2norm

        -- Sum encoding
        mat1QOp = matrixPrep matIn1norm
        mat2QOp = matrixPrep matIn2norm
        sumEncQOp = sumEncoding 1 normQOp mat1QOp mat2QOp
        matQOt = evalOp $ sumEncQOp

        normFactor = (frobeniusNormStrict matIn1) + (frobeniusNormStrict matIn2)

        finalState1 = (normFactor :+ 0) .* (apply matQOt (ket [0,0,0,0,0]))
        finalState2 = (normFactor :+ 0) .* (apply matQOt (ket [0,0,0,0,1]))
        finalState3 = (normFactor :+ 0) .* (apply matQOt (ket [0,0,0,1,0]))
     
        m11 = inner  (ket [0,0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0,0]) finalState3

        m21 = inner  (ket [0,0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,0,1]) finalState3
 
        m31 = inner  (ket [0,0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,0,1,0]) finalState3
 
        matOutData = [[m11,m12,m13], [m21,m22,m23], [m31,m32,m33]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        toLists (subMat matOut (addMat matIn1 matIn2)) --Forskellen sendes videre 


productTester :: () -> StateT--[[ComplexT]]
productTester () =
    let
        -- Matrix data
        matData1 =   [[1,2], [3,4]]
        matData2 =   [[0,1], [1,0]]
        matIn1 = fromRows matData1
        matIn2 = fromRows matData2
        
        -- Product encoding
        mat1QOp = matrixPrep matIn1
        mat2QOp = matrixPrep matIn2
        prodEncQOp = productEncoding 1 mat1QOp mat2QOp
        matQOt = evalOp $ prodEncQOp

        normFactor = (frobeniusNormStrict matIn1 :+ 0) * (frobeniusNormStrict matIn2 :+ 0)

        finalState1 = normFactor .* (apply matQOt (ket [0,0,0]))
        finalState2 = normFactor .* (apply matQOt (ket [0,0,1]))

        m11 = inner (ket [0,0,0]) finalState1
        m12 = inner (ket [0,0,0]) finalState2
        
        m21 = inner (ket [1,0,0]) finalState1 -- Jeg kan ikke få disse til at give det rigtige med 001. 
        m22 = inner (ket [1,0,0]) finalState2 -- Trods at finalState2 = ...+ 3*001 + ... får m22 IKKE værdien 3 som den bør
                                              -- Tror altså at det er en fejl som IKKE er i productEncoding

        matOutData =   [[m11,m12], [m21,m22]]

        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINGSFEJL FJERNES
    
        matInProduct = fromRows ([[2,1], [4,3]])
    in 
        finalState2--toLists (subMat matOut matInProduct) --Forskellen sendes videre 

buildRowQOpTester :: () -> [StateT] 
buildRowQOpTester =
    let 
        -- 2 dim
        resQOp2 = buildRowQOpDS (V.fromList [1:+1,2])      
        rowQOt2 = evalOp resQOp2
        -- Calculate the normalization factor and apply the operator
        finalState2 = sqrt(1^2 + 1^2 + 2^2) .* (apply rowQOt2 (ket ([0])))

        -- 3 dim
        resQOp3 = buildRowQOpDS (V.fromList [1,2,3])      
        rowQOt3 = evalOp resQOp3
        -- Calculate the normalization factor and apply the operator
        finalState3 = sqrt(1^2 + 2^2 + 3^2) .* (apply rowQOt3 (ket ([0,0])))

        -- 4 dim
        resQOp4 = buildRowQOpDS (V.fromList [1,2,3,4])      
        rowQOt4 = evalOp resQOp4
        -- Calculate the normalization factor and apply the operator
        finalState4 = sqrt(1^2 + 2^2 + 3^2 + 4^2) .* (apply rowQOt4 (ket ([0,0])))

        -- 5 dim
        resQOp5 = buildRowQOpDS (V.fromList [1,2,3,4,5])      
        rowQOt5 = evalOp resQOp5
        -- Calculate the normalization factor and apply the operator
        finalState5 = sqrt(1^2 + 2^2 + 3^2 + 4^2 + 5^2) .* (apply rowQOt5 (ket ([0,0,0])))

    in
        return [finalState2,finalState3,finalState4,finalState5]

matrixPrepTester3 :: () -> [[ComplexT]] 
matrixPrepTester3 () =
     let
        matData =  [[ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0  ], 
                    [ 3.0 :+ 0.0, 4.0 :+ 0.0 , 5.0 :+ 6.0  ],
                    [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 9.0 :+ 0.0  ]]

        matIn = fromRows matData
        matQOt = evalOp $ matrixPrep matIn

        finalState1 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,0]))
        finalState2 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,1]))
        finalState3 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,1,0]))
     
        m11 = inner  (ket [0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0]) finalState3

        m21 = inner  (ket [0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,1]) finalState3
 
        m31 = inner  (ket [0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,1,0]) finalState3
 
        matOutData = [[m11,m12,m13], [m21,m22,m23], [m31,m32,m33]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        toLists (subMat matOut matIn) --Forskellen sendes videre  



matrixPrepTester4 :: () -> [[ComplexT]] 
matrixPrepTester4 () =
     let
        matData =  [[ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0 ,  6.0 :+ 0.0 ], 
                    [ 3.0 :+ 0.0, 4.0 :+ 0.0 , 5.0 :+ 6.0 ,  7.0 :+ 0.0 ],
                    [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 9.0 :+ 0.0 ,  0.0 :+ 0.0 ],
                    [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 0.0 :+ 0.0 , 10.0 :+ 0.0 ]]

        matIn = fromRows matData
        matQOt = evalOp $ matrixPrep matIn

        finalState1 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,0]))
        finalState2 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,1]))
        finalState3 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,1,0]))
        finalState4 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,1,1]))
     
        m11 = inner  (ket [0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0]) finalState3
        m14 = inner  (ket [0,0,0,0]) finalState4

        m21 = inner  (ket [0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,1]) finalState3
        m24 = inner  (ket [0,0,0,1]) finalState4

        m31 = inner  (ket [0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,1,0]) finalState3
        m34 = inner  (ket [0,0,1,0]) finalState4
     
        m41 = inner  (ket [0,0,1,1]) finalState1
        m42 = inner  (ket [0,0,1,1]) finalState2
        m43 = inner  (ket [0,0,1,1]) finalState3
        m44 = inner  (ket [0,0,1,1]) finalState4

        matOutData = [[m11,m12,m13,m14], [m21,m22,m23,m24], [m31,m32,m33,m34],[m41,m42,m43,m44]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        toLists (subMat matOut matIn) --Forskellen sendes videre  

matrixPrepTester5 :: () -> [[ComplexT]]
matrixPrepTester5 () = 
    let
        matData =  [[ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0 ,  6.0 :+ 0.0 ,    0.0 :+ 7.0 ], 
                    [ 3.0 :+ 0.0, 4.0 :+ 0.0 , 5.0 :+ 6.0 ,  7.0 :+ 0.0 , (-8.0) :+ 0.0 ],
                    [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 9.0 :+ 0.0 ,  0.0 :+ 0.0 ,    0.0 :+ 0.0 ],
                    [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 0.0 :+ 0.0 , 10.0 :+ 0.0 ,    0.0 :+ 0.0 ],
                    [ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0 ,  0.0 :+ 6.0 ,   11.0 :+ 0.0 ]]

        matIn = fromRows matData
        matQOt = evalOp $ matrixPrep matIn

        finalState1 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,0,0,0]))
        finalState2 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,0,0,1]))
        finalState3 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,0,1,0]))
        finalState4 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,0,1,1]))
        finalState5 = (frobeniusNormStrict matIn :+ 0) .* (apply matQOt (ket [0,0,0,1,0,0]))

        m11 = inner  (ket [0,0,0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0,0,0]) finalState3
        m14 = inner  (ket [0,0,0,0,0,0]) finalState4
        m15 = inner  (ket [0,0,0,0,0,0]) finalState5

        m21 = inner  (ket [0,0,0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,0,0,1]) finalState3
        m24 = inner  (ket [0,0,0,0,0,1]) finalState4
        m25 = inner  (ket [0,0,0,0,0,1]) finalState5

        m31 = inner  (ket [0,0,0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,0,0,1,0]) finalState3
        m34 = inner  (ket [0,0,0,0,1,0]) finalState4
        m35 = inner  (ket [0,0,0,0,1,0]) finalState5
        
        m41 = inner  (ket [0,0,0,0,1,1]) finalState1
        m42 = inner  (ket [0,0,0,0,1,1]) finalState2
        m43 = inner  (ket [0,0,0,0,1,1]) finalState3
        m44 = inner  (ket [0,0,0,0,1,1]) finalState4
        m45 = inner  (ket [0,0,0,0,1,1]) finalState5

        m51 = inner  (ket [0,0,0,1,0,0]) finalState1
        m52 = inner  (ket [0,0,0,1,0,0]) finalState2
        m53 = inner  (ket [0,0,0,1,0,0]) finalState3
        m54 = inner  (ket [0,0,0,1,0,0]) finalState4
        m55 = inner  (ket [0,0,0,1,0,0]) finalState5

        matOutData = [[m11,m12,m13,m14,m15], [m21,m22,m23,m24,m25], [m31,m32,m33,m34,m35],[m41,m42,m43,m44,m45],[m51,m52,m53,m54,m55]]
        matOut = fromRows (map (map (roundComplex 9)) matOutData) -- AFRUNDINNGSFEJL
    in 
        toLists (subMat matOut matIn) --Forskellen sendes videre 


-- Helper to round a single number to n decimal places
roundTo :: RealFloat a => Int -> a -> a
roundTo n x = fromIntegral (round (x * 10^n)) / (10^n)

-- Round both parts of a complex number
roundComplex :: RealFloat a => Int -> Complex a -> Complex a
roundComplex n (r :+ i) = (roundTo n r) :+ (roundTo n i)

