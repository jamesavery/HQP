module Main where 
import HQP
import Programs.Qubitization
import Data.Complex (Complex((:+)),magnitude)
import Polynomial.Roots (roots) -- From dsp package
import Polynomial.Basic (polyadd, polysub, polymult)

main :: IO ()
main = do
    --print (show (createRealBCTest())) 
    print (show (qspVectorPhiTest()))
    --print (show (qspPolysTest()))


-----------------------------------------------------------------  
--
--   HUSK IMPLEMENTER ALLE BETINGELSER (Fra sætningerne)
--
-----------------------------------------------------------------

createRealBCTest :: () -> PolyPair
createRealBCTest () =
    -- 1. Define polynomial coefficients:
    -- The 'dsp' package expects [a0, a1, a2...]
    let 
        -- Given a real even poly A(x) with deg(A) <= 2k and A(x)>= 0 on [-1,1] 
        -- return real polys B, C as in Lemma 6
        pp = [1,0,0,0,-0.25,0, -0.5,0, -0.25]
    in
        createRealPolyPair pp

qspVectorPhiTest :: () -> [Double]
qspVectorPhiTest () =
    let 
        p = [0,0,0.5,0,0.5]
        (pC,qC) = qspPolys p
    in
        qspVectorPhi pC qC

qspPolysTest :: () -> ([ComplexT],[ComplexT])
qspPolysTest () = 
    let 
        p = [0,0,0.5,0,0.5]  -- polynomiet defineres
        (pC,qC) = qspPolys p
    in
        (pC,qC)