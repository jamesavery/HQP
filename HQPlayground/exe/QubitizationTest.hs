module Main where 
import HQP
import Programs.Qubitization
import Data.Complex (Complex((:+)),magnitude)
import Polynomial.Roots (roots) -- From dsp package
import Polynomial.Basic (polyadd, polysub, polymult)

main :: IO ()
main = do
    print (show (createRealBCTest()))
    print (show (signalVectorPhiTest()))
    print (show (qspPolysTest()))


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
        
        p = [2.0, 0, 1.0 , 0 , 4.0 ]
        q = [0,0,2]
        r = [0.0625, 0.0,(-0.5), 0.0, 1.0] -- x^4 + 0.06250000000 - 0.5000000000*x^2
        s = [0.9, 0 , -0.1] -- -0.1x^2 + 0.9
        t = [0.5 , 0 , 1] -- x^2 + 0.5
        pqr = polymult p (polymult q r)
        qr = polymult q r
        pq = polymult p q
        pr = polymult p r
        
        -- Given a real even poly A(x) with deg(A) <= 2k and A(x)>= 0 on [-1,1] 
        -- return real polys B, C as in Lemma 6
        pp = [1, 0, 0, 0, -1 ]
    in
        createRealPolyPair pp

qspVectorPhiTest :: () -> [Double]
qspVectorPhiTest () =
    let 
        p = [0.0 :+ (-1.0), 0, 1.0 :+ 1.0,0, 0, 0]
        q = [0, 0.0 :+ sqrt(2.0),0]
    in
        qspVectorPhi p q

qspPolysTest :: () -> ([ComplexT],[ComplexT])
qspPolysTest () = 
    let 
        p = [0,0,1] -- x^2 ... som skal blive til pC = (1 + i)x^2 -i (even), qC = sqrt(2)ix (odd)
                    -- Hvilket er [0.0 :+ (-1.0),0.0 :+ 0.0,1.0 :+ 1.0]           ,[0.0 :+ 0.0,0.0 :+ sqrt(2)]
                    -- Jeg får:   [0.0 :+ (-1.0),0.0 :+ 0.0,1.0 :+ 1.0,0.0 :+ 0.0],[0.0 :+ 0.0,0.0 :+ 1.4142135623730951,0.0 :+ 0.0]
        (pC,qC) = qspPolys p
    in
        (pC,qC)