module Programs.MatrixArithmetic where
import HQP
import HQP.QOp.MPSSemantics
import System.Random(mkStdGen, randoms)
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.List
import Data.Maybe
import Debug.Trace
import Data.Complex
import Data.Ratio

-- Kun Testet på reelle matricer
productEncoding :: Int -> QOp -> QOp -> QOp
productEncoding numQbits qOp1 qOp2 =
    let 
        swap = if numQbits == 1 then Permute [1,0] else Permute ([numQbits .. (2*numQbits - 1)] ++ [0 .. (numQbits - 1)])
    in
        (((Id numQbits) ⊗ (qOp1 <> swap)) ∘ (qOp2 ⊗ (Id numQbits)) ∘ ((Id numQbits) ⊗ swap))

-- Kun Testet på reelle matricer
sumEncoding :: Int -> QOp -> QOp -> QOp -> QOp
sumEncoding numQbits normQOp qOp1 qOp2 =
    let 
        wMatrix = DirectSum qOp1 qOp2
    in
        ((Adjoint normQOp) ⊗ Id (numQbits + 1)) ∘ wMatrix ∘ (normQOp ⊗ Id (numQbits + 1))