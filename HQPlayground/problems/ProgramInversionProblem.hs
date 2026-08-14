module PI where

import HQP
import Programs.QFT
import qualified HQP.QOp.MatrixSemantics as MatSem 

{-| 
   In src/Programs/QFT.hs we have implemented the Quantum Fourier Transform (QFT) as a quantum program. First have a look at the QFT implementation there.

   In this exercise, we will get comfortable working with programs as data structures that we can manipulate symbolically, performing program transformations that preserve their semantics.

   1. First generate the QFT program on 1,2, and 3 qubits, and inspect it using both show and the pretty printer. 
   
   2. Evaluate the QFT program on 1, 2, 4, 10 qubits using the matrix semantics from MatrixSemantics.hs. Verify that the resulting matrices are unitary. What would happen if you tried this with 20 qubits? Use 'cabal repl ProgramInversion' to explore this interactively.

   Look in src/HQP/QOp/Simplify.hs for an example of a program transformation that simplifies QOp syntax trees by removing Empty operators. You can use the fixEmpty function from there to simplify your QFT programs before printing or evaluating them. Check that the simplification preserves the semantics of the program.

   3. Now write a function (QOp -> QOp) that performs the inversion of a quantum program. I.e., given Adjoint (qft n), this rewrite function should actually perform the Adjoint operation recursively on all sub-operators, resulting in a new QOp that is the inverse of the original QFT program. 

   To achieve this: 1) Find out what the adjoints of primitive operators I, X, Y, X, Rz(Phi) are.  Can they be expressed as operators in our operator language?  2) What are the rules for what the adjoints of compositions, tensor products, direct sums and permutations are?  
    
   Every operator term in our operator language has an equivalent operator term without any occurrences of the  Adjoint constructor.

    4. Use this function to construct the inverse QFT program by applying it to the QFT program on n qubits. Verify that the inverse QFT program is indeed the inverse by evaluating the composition of QFT and its inverse on some test states.

    5. Write a function that takes a list of rewrite rules (e.g. [cleanEmpty,cleanAdjoint]) and applies them all to a QOp. Then write a function that applies such a list of rewrite rules repeatedly until a fixpoint is reached (i.e., applying the rules does not change the QOp anymore). You can use this to combine multiple optimization passes in your future work.
 -}

pp = putStrLn . showOp -- Pretty print QOp

prog1 = qft 1
prog2 = qft 2
prog3 = qft 3
prog4 = qft 4
prog20 = qft 20

prog3_simplified = cleanOnes prog3
inv3 = Adjoint prog3_simplified
inv4 = Adjoint (cleanOnes prog4)

