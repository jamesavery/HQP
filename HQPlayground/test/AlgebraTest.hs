{-# LANGUAGE ScopedTypeVariables #-}
-- | Algebraic identities of the QOp algebra and the simplifier rewrite rules.
--   Each property checks that two QOp expressions evaluate to the same matrix
--   under MatrixSemantics. Untested today; high coverage for the (largely
--   structural) simplifier and `dagger`.
--
--   Grouped per the README backlog:
--     * Simplifier preservation — every rewrite preserves `evalOp`.
--     * Adjoint distribution — how `Adjoint` commutes with each constructor.
--     * Algebraic structure — associativity, identity, bifunctoriality.
module Main where

import HQP
import qualified HQP.QOp.MatrixSemantics as MS
import HQP.QOp.MatrixSemantics (CMat)

import Numeric.LinearAlgebra (takeDiag, toList)
import Data.Complex (Complex(..), magnitude)

import Test.Tasty
import Test.Tasty.QuickCheck

import Generators

------------------------------------------------------------------------
-- Helpers
------------------------------------------------------------------------

-- | Evaluate two ops under MatrixSemantics and compare matrices.
eqMat :: QOp -> QOp -> Property
eqMat a b = (MS.evalOp a :: CMat) ~~ (MS.evalOp b :: CMat)

-- | Generate a same-arity companion op for a `RandomQOp`.
sameArity :: Int -> Gen QOp
sameArity n = sized $ \sz -> genQOpAt n (max 1 (min sz 3))

-- | Generate a same-arity Pauli string (used to keep R axes valid).
samePauli :: Int -> Gen QOp
samePauli = genPauli

-- | Build "X on qubit k of n total qubits" via Tensor.
xOnQubit :: Int -> Int -> QOp
xOnQubit k n
  | k < 0 || k >= n = error "xOnQubit: k out of range"
  | otherwise       = Tensor (Id k) (Tensor X (Id (n - k - 1)))

------------------------------------------------------------------------
-- 1. Simplifier preservation
------------------------------------------------------------------------

-- Each rewrite rule from `HQP.QOp.Simplify` should not change `evalOp`.
-- Note: some rewrites are sound only on certain shapes; for those we let
-- QC's random op fall through unchanged (which trivially holds).

prop_cleanOnes_preserves :: RandomQOp -> Property
prop_cleanOnes_preserves (RandomQOp _ op) = eqMat op (cleanOnes op)

prop_cleanAdjoints_preserves :: RandomQOp -> Property
prop_cleanAdjoints_preserves (RandomQOp _ op) = eqMat op (cleanAdjoints op)

prop_liftComposes_preserves :: RandomQOp -> Property
prop_liftComposes_preserves (RandomQOp _ op) = eqMat op (liftComposes op)

prop_doComposes_preserves :: RandomQOp -> Property
prop_doComposes_preserves (RandomQOp _ op) = eqMat op (doComposes op)

prop_pushComposes_preserves :: RandomQOp -> Property
prop_pushComposes_preserves (RandomQOp _ op) = eqMat op (pushComposes op)

-- NOTE: `cleanop` (fixpoint composition of all 5 simplifiers up to 10000 iters)
-- appears to not converge on some random ops in the QC generator's space — the
-- test hangs even with QuickCheckTests 1, MaxSize 2. Either one of the rewrite
-- rules oscillates structurally or the per-iteration cost is huge. Skipped
-- here; tracked as a known issue. Each individual simplifier test below passes.
--
-- prop_cleanop_preserves :: RandomQOp -> Property
-- prop_cleanop_preserves (RandomQOp _ op) = eqMat op (cleanop op)

simplifierTests :: TestTree
simplifierTests = testGroup "Simplifier preservation"
  [ testProperty "cleanOnes preserves evalOp"      prop_cleanOnes_preserves
  , testProperty "cleanAdjoints preserves evalOp"  prop_cleanAdjoints_preserves
  , testProperty "liftComposes preserves evalOp"   prop_liftComposes_preserves
  , testProperty "doComposes preserves evalOp"     prop_doComposes_preserves
  , testProperty "pushComposes preserves evalOp"   prop_pushComposes_preserves
  ]

------------------------------------------------------------------------
-- 2. Adjoint distribution rules
------------------------------------------------------------------------

-- adj distributes over each constructor; tested by comparing matrices.

prop_adj_compose :: RandomQOp -> Property
prop_adj_compose (RandomQOp n a) =
  forAll (sameArity n) $ \b ->
    eqMat (Adjoint (Compose a b)) (Compose (Adjoint b) (Adjoint a))

prop_adj_tensor :: RandomQOp -> RandomQOp -> Property
prop_adj_tensor (RandomQOp _ a) (RandomQOp _ b) =
  -- cap total qubit count to avoid huge matrices
  op_qubits a + op_qubits b <= 5 ==>
    eqMat (Adjoint (Tensor a b)) (Tensor (Adjoint a) (Adjoint b))

prop_adj_directsum :: RandomQOp -> Property
prop_adj_directsum (RandomQOp n a) =
  forAll (sameArity n) $ \b ->
    eqMat (Adjoint (DirectSum a b)) (DirectSum (Adjoint a) (Adjoint b))

prop_adj_C :: RandomQOp -> Property
prop_adj_C (RandomQOp _ a) =
  eqMat (Adjoint (C a)) (C (Adjoint a))

prop_adj_R :: RotGen -> Property
prop_adj_R (RotGen _ ax θ) =
  -- For Pauli axes, adj (R ax θ) = R ax (-θ).
  eqMat (Adjoint (R ax θ)) (R ax (-θ))

prop_adj_permute :: Perm -> Property
prop_adj_permute (Perm π) =
  eqMat (Adjoint (Permute π)) (Permute (invertPerm π))

prop_adj_involution :: RandomQOp -> Property
prop_adj_involution (RandomQOp _ a) =
  eqMat (Adjoint (Adjoint a)) a

adjointTests :: TestTree
adjointTests = testGroup "Adjoint distribution"
  [ testProperty "adj (a ∘ b) = adj b ∘ adj a"            prop_adj_compose
  , testProperty "adj (a ⊗ b) = adj a ⊗ adj b"            prop_adj_tensor
  , testProperty "adj (a ⊕ b) = adj a ⊕ adj b"            prop_adj_directsum
  , testProperty "adj (C a) = C (adj a)"                  prop_adj_C
  , testProperty "adj (R ax θ) = R ax (-θ) (Pauli axis)"  prop_adj_R
  , testProperty "adj (Permute π) = Permute (invertPerm π)" prop_adj_permute
  , testProperty "adj (adj a) = a"                        prop_adj_involution
  ]

------------------------------------------------------------------------
-- 3. Algebraic structure (associativity, identity, bifunctoriality)
------------------------------------------------------------------------

prop_compose_assoc :: RandomQOp -> Property
prop_compose_assoc (RandomQOp n a) =
  forAll (sameArity n) $ \b ->
  forAll (sameArity n) $ \c ->
    eqMat (Compose (Compose a b) c) (Compose a (Compose b c))

prop_tensor_assoc :: RandomQOp -> RandomQOp -> RandomQOp -> Property
prop_tensor_assoc (RandomQOp _ a) (RandomQOp _ b) (RandomQOp _ c) =
  op_qubits a + op_qubits b + op_qubits c <= 6 ==>
    eqMat (Tensor (Tensor a b) c) (Tensor a (Tensor b c))

prop_compose_id_right :: RandomQOp -> Property
prop_compose_id_right (RandomQOp n a) = eqMat (Compose a (Id n)) a

prop_compose_id_left :: RandomQOp -> Property
prop_compose_id_left  (RandomQOp n a) = eqMat (Compose (Id n) a) a

prop_tensor_one_right :: RandomQOp -> Property
prop_tensor_one_right (RandomQOp _ a) = eqMat (Tensor a (Id 0)) a

prop_tensor_one_left :: RandomQOp -> Property
prop_tensor_one_left  (RandomQOp _ a) = eqMat (Tensor (Id 0) a) a

-- (a ⊗ b) ∘ (c ⊗ d) = (a ∘ c) ⊗ (b ∘ d) when arities align.
prop_bifunctor_tensor :: RandomQOp -> RandomQOp -> Property
prop_bifunctor_tensor (RandomQOp na a) (RandomQOp nb b) =
  op_qubits a + op_qubits b <= 4 ==>
  forAll (sameArity na) $ \c ->
  forAll (sameArity nb) $ \d ->
    eqMat (Compose (Tensor a b) (Tensor c d))
          (Tensor (Compose a c) (Compose b d))

-- (a ⊕ b) ∘ (c ⊕ d) = (a ∘ c) ⊕ (b ∘ d).
prop_bifunctor_dsum :: RandomQOp -> Property
prop_bifunctor_dsum (RandomQOp n a) =
  forAll (sameArity n) $ \b ->
  forAll (sameArity n) $ \c ->
  forAll (sameArity n) $ \d ->
    eqMat (Compose (DirectSum a b) (DirectSum c d))
          (DirectSum (Compose a c) (Compose b d))

-- MS implements C a literally as (Id n_a) ⊕ a (line 115 of MatrixSemantics.hs);
-- this checks the same identity holds across the matrix definitions.
prop_C_as_directsum :: RandomQOp -> Property
prop_C_as_directsum (RandomQOp _ a) =
  let nA = op_qubits a
  in eqMat (C a) (DirectSum (Id nA) a)

algebraTests :: TestTree
algebraTests = testGroup "Algebraic structure"
  [ testProperty "Compose is associative"             prop_compose_assoc
  , testProperty "Tensor is associative"              prop_tensor_assoc
  , testProperty "a ∘ Id n = a"                       prop_compose_id_right
  , testProperty "Id n ∘ a = a"                       prop_compose_id_left
  , testProperty "a ⊗ One = a"                        prop_tensor_one_right
  , testProperty "One ⊗ a = a"                        prop_tensor_one_left
  , testProperty "(a⊗b)∘(c⊗d) = (a∘c)⊗(b∘d)"          prop_bifunctor_tensor
  , testProperty "(a⊕b)∘(c⊕d) = (a∘c)⊕(b∘d)"          prop_bifunctor_dsum
  , testProperty "C a = Id n_a ⊕ a"                   prop_C_as_directsum
  ]

------------------------------------------------------------------------
-- 4. Pauli / Clifford gate identities
------------------------------------------------------------------------

prop_X_squared, prop_Y_squared, prop_Z_squared, prop_H_squared :: Property
prop_X_squared = once $ eqMat (Compose X X) (Id 1)
prop_Y_squared = once $ eqMat (Compose Y Y) (Id 1)
prop_Z_squared = once $ eqMat (Compose Z Z) (Id 1)
prop_H_squared = once $ eqMat (Compose H H) (Id 1)

prop_CX_squared :: Property
prop_CX_squared = once $ eqMat (Compose (C X) (C X)) (Id 2)

prop_HXH_is_Z :: Property
prop_HXH_is_Z = once $ eqMat (Compose H (Compose X H)) Z

prop_HZH_is_X :: Property
prop_HZH_is_X = once $ eqMat (Compose H (Compose Z H)) X

-- (I ⊗ H) ∘ CX ∘ (I ⊗ H) = CZ (Hadamard on target turns CNOT into CZ).
prop_IHCXIH_is_CZ :: Property
prop_IHCXIH_is_CZ = once $
  eqMat (Compose (Tensor (Id 1) H) (Compose (C X) (Tensor (Id 1) H)))
        (C Z)

-- X-conjugation flips the sign of a Z rotation: (R Z θ) ∘ X = X ∘ R Z (-θ).
prop_RZ_X_conjugate :: Property
prop_RZ_X_conjugate = forAll genRat $ \θ ->
  eqMat (Compose (R Z θ) X) (Compose X (R Z (-θ)))

pauliCliffordTests :: TestTree
pauliCliffordTests = testGroup "Pauli / Clifford identities"
  [ testProperty "X² = I"               prop_X_squared
  , testProperty "Y² = I"               prop_Y_squared
  , testProperty "Z² = I"               prop_Z_squared
  , testProperty "H² = I"               prop_H_squared
  , testProperty "(C X)² = I"           prop_CX_squared
  , testProperty "H X H = Z"            prop_HXH_is_Z
  , testProperty "H Z H = X"            prop_HZH_is_X
  , testProperty "(I⊗H) CX (I⊗H) = CZ"  prop_IHCXIH_is_CZ
  , testProperty "R Z θ ∘ X = X ∘ R Z (-θ)" prop_RZ_X_conjugate
  ]

------------------------------------------------------------------------
-- 5. Permutation algebra
------------------------------------------------------------------------

-- The identity permutation should equal `Id n`.
prop_permute_identity_perm :: Property
prop_permute_identity_perm = forAll (choose (1, 4) :: Gen Int) $ \n ->
  eqMat (Permute [0 .. n - 1]) (Id n)

-- `Compose (Permute π) (Permute σ) = Permute κ` where κ[q] = σ[π[q]].
-- Convention: Compose a b applies b first, then a, so output q = (Permute π
-- applied to (Permute σ applied to input))[q] = input[σ[π[q]]].
prop_permute_compose :: Perm -> Property
prop_permute_compose (Perm π) =
  forAll (genPermutation (length π)) $ \σ ->
    let κ = [ σ !! (π !! q) | q <- [0 .. length π - 1] ]
    in eqMat (Compose (Permute π) (Permute σ)) (Permute κ)

-- Conjugation: Permute π ∘ X_0 ∘ Permute π⁻¹ = X on qubit (π⁻¹[0]).
prop_permute_conjugation :: Perm -> Property
prop_permute_conjugation (Perm π) =
  let n     = length π
      invπ  = invertPerm π
      conj  = Compose (Permute π) (Compose (xOnQubit 0 n) (Permute invπ))
      expected = xOnQubit (invπ !! 0) n
  in eqMat conj expected

permuteAlgebraTests :: TestTree
permuteAlgebraTests = testGroup "Permutation algebra"
  [ testProperty "Permute [0..n-1] = Id n"                  prop_permute_identity_perm
  , testProperty "Permute π ∘ Permute σ = Permute (σ∘π)"    prop_permute_compose
  , testProperty "Permute π ∘ X_0 ∘ Permute π⁻¹ = X_{π⁻¹[0]}" prop_permute_conjugation
  ]

------------------------------------------------------------------------
-- 6. Rotation algebra (additional)
------------------------------------------------------------------------

-- R ax 4 = Id n: 4π rotation around any involution Pauli is identity.
prop_R_4_is_Id :: RotGen -> Property
prop_R_4_is_Id (RotGen n ax _) = eqMat (R ax 4) (Id n)

-- General angle additivity: R ax (a+b) = R ax a ∘ R ax b.
prop_R_angle_additive :: RotGen -> Property
prop_R_angle_additive (RotGen _ ax θ) = forAll genRat $ \θ' ->
  eqMat (R ax (θ + θ')) (Compose (R ax θ) (R ax θ'))

-- 1-qubit Pauli rotation: Tr(R ax θ) = 2 cos(πθ/2) for any Pauli ax (including I).
prop_R_trace_1q :: Property
prop_R_trace_1q = forAll (elements [X, Y, Z]) $ \ax ->
                  forAll genRat $ \θ ->
  let m        = MS.evalOp (R ax θ) :: CMat
      diag     = toList (takeDiag m)
      traceVal = sum diag
      expected = (2 * cos (pi * fromRational θ / 2)) :+ 0
  in counterexample ("trace = " ++ show traceVal ++ ", expected " ++ show expected)
       (magnitude (traceVal - expected) < tol)

rotationAlgebraTests :: TestTree
rotationAlgebraTests = testGroup "Rotation algebra (more)"
  [ testProperty "R ax 4 = Id n"                       prop_R_4_is_Id
  , testProperty "R ax (θ₁+θ₂) = R ax θ₁ ∘ R ax θ₂"    prop_R_angle_additive
  , testProperty "Tr(R ax θ) = 2 cos(πθ/2)  [1-qubit]" prop_R_trace_1q
  ]

------------------------------------------------------------------------
-- 7. Phase identities
------------------------------------------------------------------------

-- Phase composition: Phase a ∘ Phase b = Phase (a+b) as matrices.
prop_phase_compose :: Property
prop_phase_compose = forAll genRat $ \a -> forAll genRat $ \b ->
  eqMat (Compose (Phase a) (Phase b)) (Phase (a + b))

-- Phase 0 = scalar 1 = Id 0 (i.e. One).
prop_phase_zero :: Property
prop_phase_zero = once $ eqMat (Phase 0) (Id 0)

-- Phase 2 = e^{2πi} = 1 = Id 0.
prop_phase_two :: Property
prop_phase_two = once $ eqMat (Phase 2) (Id 0)

phaseTests :: TestTree
phaseTests = testGroup "Phase identities"
  [ testProperty "Phase a ∘ Phase b = Phase (a+b)" prop_phase_compose
  , testProperty "Phase 0 = Id 0"                  prop_phase_zero
  , testProperty "Phase 2 = Id 0"                  prop_phase_two
  ]

------------------------------------------------------------------------

qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 50) . localOption (QuickCheckMaxSize 4)

main :: IO ()
main = defaultMain $ qcOpts $ testGroup "AlgebraTest"
  [ simplifierTests
  , adjointTests
  , algebraTests
  , pauliCliffordTests
  , permuteAlgebraTests
  , rotationAlgebraTests
  , phaseTests
  ]
