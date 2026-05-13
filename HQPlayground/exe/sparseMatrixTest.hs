module Main where

import HQP.QOp
import HQP.QOp.StatevectorSemantics
import HQP.PrettyPrint
import Programs.QFT

import Data.Ratio


main :: IO ()
main = do
    --printS $ (permutationTest ()) 
    --printS $ (qftPermTest ())
    printS $ (phaseTest ())
    

permutationTest :: () -> StateT
permutationTest () =
    let 
        a = Permute [1,0,3,2]
    in
        apply (evalOp a)  (1 .* ket [1,0,0,0] .+ 2 .* ket [0,1,0,0] .+ 3 .* ket [0,0,1,0] .+ 4 .* ket [0,0,0,1])

phaseTest :: () -> StateT
phaseTest () =
    let 
        p = Phase (1 % 2)
        --q = R Z (1 % 4)
        phi0 = apply (evalOp (C p)) (ket [0,0])
        --phi1 = apply (evalOp q) phi0
    in
        phi0

qftPermTest :: () -> StateT
qftPermTest () =
    let 
        p = Permute [2,1,0]
        phi0 = apply (evalOp p) (ket [0,0,0])
        
        qftQOp = myQFT3 ()
        qftQOt = evalOp (qftQOp)
        phi1 = apply qftQOt phi0

        rotQOp = (R Z (1 % 4)) ⊗ (R Z (1 % 2)) ⊗ (R Z 1)
        rotQOt = evalOp rotQOp
        phi2 = apply rotQOt phi1

        adjQftQOt = evalOp (Adjoint (myQFT3 ()))
        phi3 = apply adjQftQOt phi2
    in
        apply (evalOp p) phi3



myQFT3 :: () -> QOp
myQFT3 () = 
    let 
        rotK = R Z (1 % 2) -- <> Phase (1 % 4)) 

        --layer1 = (H ⊗ I ⊗ I)
        layer2 = (C rotK) ⊗ I 
        --layer3 = Permute [0,2,1]
        --layer4 = (C (Phase (1%8) <> R Z (1%4)) ⊗ I)
        --layer5 = Permute [0,2,1]
        --layer6 = (I ⊗ H ⊗ I)
        --layer7 = (I ⊗ C (Phase (1%4) <> R Z (1%2)))
        --layer8 = (I ⊗ C (Phase (1%8) <> R Z (1%4)))
        --layer9 = Permute [2,1,0]
    in
        --layer9 <> layer8 <> layer7 <> layer6 <> layer5 <> layer4 <> layer3 <> layer2 <> layer1
        --layer9 <> layer6 <> layer5 <> layer3 <> layer2 <> layer1
        layer2


{- Output fra qft 2

Compose 
(
    Permute [1,0]
) 
(
    Compose 
    (
        Tensor H (Id 1)
    ) 
    (
        Compose 
        (
            Compose 
            (
                Tensor 
                (
                    Id 0
                ) 
                (
                    C 
                    (
                        Tensor 
                        (
                            Id 0
                        ) 
                        (
                            R Z (1 % 2)
                        )
                    )
                )
            ) 
            (
                Id 2
            )
        ) 
        (
            Tensor 
            (
                Id 1
            ) 
                H
        )
    )
)
-}
