module Main where

import HQP.QOp
import HQP.QOp.StatevectorSemantics
import HQP.PrettyPrint
import Programs.QFT
import Programs.SparseMatrixPreparation
import Data.Ratio
import qualified Data.Vector as V
import Numeric.LinearAlgebra (fromLists, toLists, norm_Frob)
import Data.Complex


main :: IO ()
main = do
    --print $ (show (permutationTest ())) 
    --printS $ (qftPermTest ())
    print $ (show (sparseUAEncodingTest4 ()))
    print $ (show (sparseUAEncodingTest8 ()))

sparseUAEncodingTest4 :: () -> [[ComplexT]]
sparseUAEncodingTest4 () =
    let
        sparseValues = [[1], [-1]]
        valOracleQOp = valueOracle $ buildRotVec 4 sparseValues

        permutations = V.fromList [ Permute [1,0], X⊕X  ]
        posOracleQOp = positionOracle permutations

        uaQOp = sparseUAEncoding 1 valOracleQOp posOracleQOp
        uaQOt = evalOp uaQOp 

        finalState1 =  (apply uaQOt (ket [0,0,0,0]))
        finalState2 =  (apply uaQOt (ket [0,0,0,1]))
        finalState3 =  (apply uaQOt (ket [0,0,1,0]))
        finalState4 =  (apply uaQOt (ket [0,0,1,1]))

        m11 = inner  (ket [0,0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0,0]) finalState3
        m14 = inner  (ket [0,0,0,0,0]) finalState4

        m21 = inner  (ket [0,0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,0,1]) finalState3
        m24 = inner  (ket [0,0,0,0,1]) finalState4

        m31 = inner  (ket [0,0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,0,1,0]) finalState3
        m34 = inner  (ket [0,0,0,1,0]) finalState4

        m41 = inner  (ket [0,0,0,1,1]) finalState1
        m42 = inner  (ket [0,0,0,1,1]) finalState2
        m43 = inner  (ket [0,0,0,1,1]) finalState3
        m44 = inner  (ket [0,0,0,1,1]) finalState4

        matOutData =   [[m11,m12,m13,m14], 
                        [m21,m22,m23,m24], 
                        [m31,m32,m33,m34],
                        [m41,m42,m43,m44]
                        ]
        matOut = fromLists matOutData
    in 
        matOutData
        
sparseUAEncodingTest8 :: () -> [[ComplexT]]
sparseUAEncodingTest8 () =
    let
        sparseValues = [[1], [-1]]
        valOracleQOp = valueOracle $ buildRotVec 8 sparseValues

        permutations = V.fromList [ Permute [1,2,0],  (( (X ⊕ I) ⊕ (Id 2)) )<>(( Id 2)⊕ (I ⊕ X)) ]
        posOracleQOp = positionOracle permutations

        uaQOp = sparseUAEncoding 1 valOracleQOp posOracleQOp
        uaQOt = evalOp uaQOp 

        finalState1 =  (apply uaQOt (ket [0,0,0,0,0]))
        finalState2 =  (apply uaQOt (ket [0,0,0,0,1]))
        finalState3 =  (apply uaQOt (ket [0,0,0,1,0]))
        finalState4 =  (apply uaQOt (ket [0,0,0,1,1]))
        finalState5 =  (apply uaQOt (ket [0,0,1,0,0]))
        finalState6 =  (apply uaQOt (ket [0,0,1,0,1]))
        finalState7 =  (apply uaQOt (ket [0,0,1,1,0]))
        finalState8 =  (apply uaQOt (ket [0,0,1,1,1]))

        m11 = inner  (ket [0,0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0,0]) finalState3
        m14 = inner  (ket [0,0,0,0,0]) finalState4
        m15 = inner  (ket [0,0,0,0,0]) finalState5
        m16 = inner  (ket [0,0,0,0,0]) finalState6
        m17 = inner  (ket [0,0,0,0,0]) finalState7
        m18 = inner  (ket [0,0,0,0,0]) finalState8

        m21 = inner  (ket [0,0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,0,1]) finalState3
        m24 = inner  (ket [0,0,0,0,1]) finalState4
        m25 = inner  (ket [0,0,0,0,1]) finalState5
        m26 = inner  (ket [0,0,0,0,1]) finalState6
        m27 = inner  (ket [0,0,0,0,1]) finalState7
        m28 = inner  (ket [0,0,0,0,1]) finalState8

        m31 = inner  (ket [0,0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,0,1,0]) finalState3
        m34 = inner  (ket [0,0,0,1,0]) finalState4
        m35 = inner  (ket [0,0,0,1,0]) finalState5
        m36 = inner  (ket [0,0,0,1,0]) finalState6
        m37 = inner  (ket [0,0,0,1,0]) finalState7
        m38 = inner  (ket [0,0,0,1,0]) finalState8

        m41 = inner  (ket [0,0,0,1,1]) finalState1
        m42 = inner  (ket [0,0,0,1,1]) finalState2
        m43 = inner  (ket [0,0,0,1,1]) finalState3
        m44 = inner  (ket [0,0,0,1,1]) finalState4
        m45 = inner  (ket [0,0,0,1,1]) finalState5
        m46 = inner  (ket [0,0,0,1,1]) finalState6
        m47 = inner  (ket [0,0,0,1,1]) finalState7
        m48 = inner  (ket [0,0,0,1,1]) finalState8

        m51 = inner  (ket [0,0,1,0,0]) finalState1
        m52 = inner  (ket [0,0,1,0,0]) finalState2
        m53 = inner  (ket [0,0,1,0,0]) finalState3
        m54 = inner  (ket [0,0,1,0,0]) finalState4
        m55 = inner  (ket [0,0,1,0,0]) finalState5
        m56 = inner  (ket [0,0,1,0,0]) finalState6
        m57 = inner  (ket [0,0,1,0,0]) finalState7
        m58 = inner  (ket [0,0,1,0,0]) finalState8

        m61 = inner  (ket [0,0,1,0,1]) finalState1
        m62 = inner  (ket [0,0,1,0,1]) finalState2
        m63 = inner  (ket [0,0,1,0,1]) finalState3
        m64 = inner  (ket [0,0,1,0,1]) finalState4
        m65 = inner  (ket [0,0,1,0,1]) finalState5
        m66 = inner  (ket [0,0,1,0,1]) finalState6
        m67 = inner  (ket [0,0,1,0,1]) finalState7
        m68 = inner  (ket [0,0,1,0,1]) finalState8

        m71 = inner  (ket [0,0,1,1,0]) finalState1
        m72 = inner  (ket [0,0,1,1,0]) finalState2
        m73 = inner  (ket [0,0,1,1,0]) finalState3
        m74 = inner  (ket [0,0,1,1,0]) finalState4
        m75 = inner  (ket [0,0,1,1,0]) finalState5
        m76 = inner  (ket [0,0,1,1,0]) finalState6
        m77 = inner  (ket [0,0,1,1,0]) finalState7
        m78 = inner  (ket [0,0,1,1,0]) finalState8

        m81 = inner  (ket [0,0,1,1,1]) finalState1
        m82 = inner  (ket [0,0,1,1,1]) finalState2
        m83 = inner  (ket [0,0,1,1,1]) finalState3
        m84 = inner  (ket [0,0,1,1,1]) finalState4
        m85 = inner  (ket [0,0,1,1,1]) finalState5
        m86 = inner  (ket [0,0,1,1,1]) finalState6
        m87 = inner  (ket [0,0,1,1,1]) finalState7
        m88 = inner  (ket [0,0,1,1,1]) finalState8

        matOutData =   [[m11,m12,m13,m14,m15,m16,m17,m18], 
                        [m21,m22,m23,m24,m25,m26,m27,m28], 
                        [m31,m32,m33,m34,m35,m36,m37,m38],
                        [m41,m42,m43,m44,m45,m46,m47,m48],
                        [m51,m52,m53,m54,m55,m56,m57,m58], 
                        [m61,m62,m63,m64,m65,m66,m67,m68], 
                        [m71,m72,m73,m74,m75,m76,m77,m78],
                        [m81,m82,m83,m84,m85,m86,m87,m88]
                        ]
    in 
        roundMatrix matOutData

-- | Round a nested list matrix of complex numbers to 16 decimal places
roundMatrix :: [[ComplexT]] -> [[ComplexT]]
roundMatrix = map (map roundComplex)

-- | Round an individual complex number component to 16 decimal places
roundComplex :: ComplexT -> ComplexT
roundComplex (r :+ i) = roundVal r :+ roundVal i
  where
    factor = 10 ^ (12 :: Int)  -- 10^16 scaling metric
    roundVal x = fromIntegral (round (x * factor) :: Integer) / factor

permutationTest :: () -> [[ComplexT]]
permutationTest () =
    let 
        a = (( (X ⊕ I) ⊕ (Id 2)) )<>(( Id 2)⊕ (I ⊕ X))
        at = evalOp a
        
        finalState1 =  (apply at (ket [0,0,0]))
        finalState2 =  (apply at (ket [0,0,1]))
        finalState3 =  (apply at (ket [0,1,0]))
        finalState4 =  (apply at (ket [0,1,1]))
        finalState5 =  (apply at (ket [1,0,0]))
        finalState6 =  (apply at (ket [1,0,1]))
        finalState7 =  (apply at (ket [1,1,0]))
        finalState8 =  (apply at (ket [1,1,1]))

        m11 = inner  (ket [0,0,0]) finalState1
        m12 = inner  (ket [0,0,0]) finalState2
        m13 = inner  (ket [0,0,0]) finalState3
        m14 = inner  (ket [0,0,0]) finalState4
        m15 = inner  (ket [0,0,0]) finalState5
        m16 = inner  (ket [0,0,0]) finalState6
        m17 = inner  (ket [0,0,0]) finalState7
        m18 = inner  (ket [0,0,0]) finalState8

        m21 = inner  (ket [0,0,1]) finalState1
        m22 = inner  (ket [0,0,1]) finalState2
        m23 = inner  (ket [0,0,1]) finalState3
        m24 = inner  (ket [0,0,1]) finalState4
        m25 = inner  (ket [0,0,1]) finalState5
        m26 = inner  (ket [0,0,1]) finalState6
        m27 = inner  (ket [0,0,1]) finalState7
        m28 = inner  (ket [0,0,1]) finalState8

        m31 = inner  (ket [0,1,0]) finalState1
        m32 = inner  (ket [0,1,0]) finalState2
        m33 = inner  (ket [0,1,0]) finalState3
        m34 = inner  (ket [0,1,0]) finalState4
        m35 = inner  (ket [0,1,0]) finalState5
        m36 = inner  (ket [0,1,0]) finalState6
        m37 = inner  (ket [0,1,0]) finalState7
        m38 = inner  (ket [0,1,0]) finalState8

        m41 = inner  (ket [0,1,1]) finalState1
        m42 = inner  (ket [0,1,1]) finalState2
        m43 = inner  (ket [0,1,1]) finalState3
        m44 = inner  (ket [0,1,1]) finalState4
        m45 = inner  (ket [0,1,1]) finalState5
        m46 = inner  (ket [0,1,1]) finalState6
        m47 = inner  (ket [0,1,1]) finalState7
        m48 = inner  (ket [0,1,1]) finalState8

        m51 = inner  (ket [1,0,0]) finalState1
        m52 = inner  (ket [1,0,0]) finalState2
        m53 = inner  (ket [1,0,0]) finalState3
        m54 = inner  (ket [1,0,0]) finalState4
        m55 = inner  (ket [1,0,0]) finalState5
        m56 = inner  (ket [1,0,0]) finalState6
        m57 = inner  (ket [1,0,0]) finalState7
        m58 = inner  (ket [1,0,0]) finalState8

        m61 = inner  (ket [1,0,1]) finalState1
        m62 = inner  (ket [1,0,1]) finalState2
        m63 = inner  (ket [1,0,1]) finalState3
        m64 = inner  (ket [1,0,1]) finalState4
        m65 = inner  (ket [1,0,1]) finalState5
        m66 = inner  (ket [1,0,1]) finalState6
        m67 = inner  (ket [1,0,1]) finalState7
        m68 = inner  (ket [1,0,1]) finalState8

        m71 = inner  (ket [1,1,0]) finalState1
        m72 = inner  (ket [1,1,0]) finalState2
        m73 = inner  (ket [1,1,0]) finalState3
        m74 = inner  (ket [1,1,0]) finalState4
        m75 = inner  (ket [1,1,0]) finalState5
        m76 = inner  (ket [1,1,0]) finalState6
        m77 = inner  (ket [1,1,0]) finalState7
        m78 = inner  (ket [1,1,0]) finalState8

        m81 = inner  (ket [1,1,1]) finalState1
        m82 = inner  (ket [1,1,1]) finalState2
        m83 = inner  (ket [1,1,1]) finalState3
        m84 = inner  (ket [1,1,1]) finalState4
        m85 = inner  (ket [1,1,1]) finalState5
        m86 = inner  (ket [1,1,1]) finalState6
        m87 = inner  (ket [1,1,1]) finalState7
        m88 = inner  (ket [1,1,1]) finalState8

        matOutData =   [[m11,m12,m13,m14,m15,m16,m17,m18], 
                        [m21,m22,m23,m24,m25,m26,m27,m28], 
                        [m31,m32,m33,m34,m35,m36,m37,m38],
                        [m41,m42,m43,m44,m45,m46,m47,m48],
                        [m51,m52,m53,m54,m55,m56,m57,m58], 
                        [m61,m62,m63,m64,m65,m66,m67,m68], 
                        [m71,m72,m73,m74,m75,m76,m77,m78],
                        [m81,m82,m83,m84,m85,m86,m87,m88]
                        ]
    in 
        roundMatrix matOutData


qftPermTest :: () -> StateT
qftPermTest () =
    let 
        p = Permute [2,1,0]
        phi0 = apply (evalOp p) (ket [0,0,0])
        
        qftQOp = qft 3
        qftQOt = evalOp (qftQOp)
        phi1 = apply qftQOt phi0

        rotQOp = (R Z (1 % 4)) ⊗ (R Z (1 % 2)) ⊗ (R Z 1)
        rotQOt = evalOp rotQOp
        phi2 = apply rotQOt phi1

        adjQftQOt = evalOp (Adjoint (qft 3))
        phi3 = apply adjQftQOt phi2
    in
        apply (evalOp p) phi3