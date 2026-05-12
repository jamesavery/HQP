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
sumEncoding :: QOp -> QOp -> QOp -> QOp
sumEncoding normQOp qOp1 qOp2 =
    let 
        wMatrix = DirectSum qOp1 qOp2
        n = (op_qubits wMatrix) - (op_qubits normQOp)
    in
        ((Adjoint normQOp) ⊗ (Id n)) ∘ wMatrix ∘ (normQOp ⊗ (Id n))







{-


(
(
    ( 
        adj 
        (
            (R Z (0 % 1) ∘ (R Y ((-5788800705379855) % 9007199254740992) ∘ R Z (0 % 1))) ⊗ One
        )
    ) 
    ⊗ 
    Id 2
) 
∘ 
(
    (
        (
            (
                adj 
                (
                    (
                        (
                            (
                                (
                                    (
                                        (R Z (0 % 1) ∘ (R Y ((-1587142288228767) % 2251799813685248) ∘ R Z (0 % 1))) 
                                        ⊕ 
                                        (   
                                            R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1))
                                        )
                                    ) 
                                    ⊗ 
                                    One
                                ) 
                                ∘ 
                                (
                                    R Y ((-5334341100569611) % 9007199254740992) ⊗ I
                                )
                            ) 
                            ⊕ 
                            (
                                (
                                    (
                                        (R Z (0 % 1) ∘ (R Y ((-5138125964800215) % 9007199254740992) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-4317294479762133) % 9007199254740992) ⊗ I))) ⊕ (((((R Z (0 % 1) ∘ (R Y ((-2442656102601617) % 4503599627370496) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-8057015789365283) % 18014398509481984) ⊗ I)) ⊕ Id 2)) ∘ Permute [2,3,0,1])) ∘ (((((R Z (0 % 1) ∘ (R Y ((-837000630181079) % 1125899906842624) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-5563767763481133) % 9007199254740992) ⊗ I)) ⊗ Id 2)) ⊕ ((adj (((((((R Z (0 % 1) ∘ (R Y ((-4776448809064993) % 9007199254740992) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-7788369832850731) % 18014398509481984) ⊗ I)) ⊕ ((((R Z (0 % 1) ∘ (R Y ((-4715878937184305) % 9007199254740992) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-3816723671856269) % 9007199254740992) ⊗ I))) ⊕ (((((R Z (0 % 1) ∘ (R Y ((-2338654464128839) % 4503599627370496) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-1883179745058753) % 4503599627370496) ⊗ I)) ⊕ Id 2)) ∘ Permute [2,3,0,1])) ∘ (((((R Z (0 % 1) ∘ (R Y ((-2592750248214573) % 4503599627370496) ∘ R Z (0 % 1))) ⊕ (R Z (0 % 1) ∘ (R Y (0 % 1) ∘ R Z (0 % 1)))) ⊗ One) ∘ (R Y ((-8736738316232525) % 18014398509481984) ⊗ I)) ⊗ Id 2))) ∘ (((R Z (0 % 1) ∘ (R Y ((-5788800705379855) % 9007199254740992) ∘ R Z (0 % 1))) ⊗ One) ⊗ Id 2)))

-}