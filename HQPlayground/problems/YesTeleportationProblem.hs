module Main where 

import QubitOperators
import MatrixSemantics
import PrettyMatrix
import PrettyOp
import System.Random(mkStdGen, randoms)
-- import Data.List(intercalate)

{-| Exercise: Implement the quantum teleportation protocol from N&C Section 1.3.7 (Fig. 1.11), and demonstrate that teleportation of the full quantum state, contrary to cloning, *is* possible.

Use the operations in QubitsOperator.hs.

- How do you implement X^M2, Z^M1 without dynamic control? (Hint: A^0 = I for any operator A).

- Our simple language doesn't have qubit numbering. How can you implement the flow from measurement M1 to qubit 2? (Hint: op can be multi-qubit in our controlled operation 'C op')
-}
teleprog :: Program
teleprog = let
        op1 = Unitary $ 
            I ⊗ I ⊗ I -- Replace by your own program
        
        op2 = Unitary $ 
            I ⊗ I ⊗ I -- Replace by your own program
    in
        [op1, Measure [0,1], op2]

{-| 'teleport rng psi' runs the teleprog program on the input 1Q-state psi. rng is a
list of numbers [0..1] for random measurements. See below for how to make an infinite list of random numbers, or supply a fixed list for debugging. 

returns: (final_state, rng') where rng' 
is the updated random number generator/list.
-}
teleport :: RNG -> StateT -> (StateT, RNG)
teleport rng psi = let
        bell    =  ket[0,0] -- replace by actual definition, N&C Section 1.3.6
        psi_bell = psi ⊗ bell
    in
        evalProg rng teleprog psi_bell


main :: IO()
main = do
    let rng0 = randoms (mkStdGen 42) :: [Double]    
    --let rng0 = [0,0,0] -- (always measure 0)
    -- let rng0 = [0,1,0] -- (First meausure 0, then measure 1)
    
    putStr $ "\nTeleport program:\n" ++ showProgram teleprog ++ "\n\n"

--    putStr $ "|ψ>   = "++(showState psi) ++ "\n" -- ". Running teleport program!\n"
    --putStr $ "|ψbb> = "++(showState psi_bell)++".\nRunning teleport program!\n"

    -- let (end_state,_) = teleport rng0 psi

    -- putStr $ "Final 3-qubit state:\n" ++ (showState end_state) ++ "\n\n"

    

    
