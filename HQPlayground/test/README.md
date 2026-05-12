# test/

Validation test suites for HQP. Built with `tasty` + `tasty-hunit` +
`tasty-quickcheck`. Shared deps live in the `test-common` cabal stanza.

## Test exes

| Exe | What it covers |
|---|---|
| `BackendValidation` | Cross-backend matrix equivalence (MatrixSemantics ↔ MPSSemantics ↔ StatevectorSemantics), unitarity, adjoint involution, rotation algebra, permutation inverse. Has caught six MPS/SV bugs to date. |
| `MatrixPrepTest` | Pure helpers in `Programs.MatrixPreparation` (`toBitString`, `splitList`, `padToPowerOf2`); `buildRowQOp` row preparation; `matrixPrep` block encoding correctness and unitarity; three-backend agreement on prepared circuits. |

Run individually: `cabal run BackendValidation` or `cabal run MatrixPrepTest`.

## Generators (in `BackendValidation.hs`)

| Generator | Produces |
|---|---|
| `RandomQOp` | An `n`-qubit `QOp` with recursion depth ≤ `d` (both ∈ [1,4]). Arity-preserving shrink. |
| `genQOpAt n d` | Same, parameterised — used inside `RandomQOp`'s `arbitrary`. |
| `genQOpBase n` | A leaf op of arity `n`: primitive (X, Y, Z, H, SX, Id), `Phase q`, `Permute π`, or `R ax θ`. |
| `genPauli n` | An `n`-qubit Pauli string (tensor of `Id 1`/`X`/`Y`/`Z`). Used as `R` axis. |
| `genPermutation n` | A permutation of `[0..n-1]`, via QuickCheck's `shuffle`. |
| `genRat` | A rational in `[-7,7]/[1,8]`. |
| `RandomVec`, `RandomMatrix` | Random complex vectors / square matrices (small integer entries). In `MatrixPrepTest`. |

QC budget: `QuickCheckTests 50`, `QuickCheckMaxSize 5` (4 in MatrixPrepTest).
The shrinkers are arity-preserving so the shrunk counterexample is still a
valid input for the property.

## Properties asserted today

### `BackendValidation`

| Property | Idea |
|---|---|
| `MPS matches MS`, `SV matches MS` | Full-matrix equivalence via `to (MPS.evalOp op) :: CMat` vs `MS.evalOp op` (and likewise SV). |
| `MPS output is unitary` | `U · U† = I` for the materialized MPS matrix. |
| `op ∘ adj op = I (MPS)` | End-to-end inverse. |
| `adj² (op) = op (MPS)` | Adjoint involution. |
| `R ax 0 = Id` | Rotation zero point. |
| `R ax 2θ = (R ax θ)²` | Angle additivity. |
| `R ax θ ∘ R ax (-θ) = I` | Rotation inverse. |
| `Permute π ∘ π⁻¹ = Id` | Permutation inverse. |

### `MatrixPrepTest`

| Property | Idea |
|---|---|
| Block extraction | Upper-left `2^n × 2^n` of `evalOp (matrixPrep M)` equals `M / ‖M‖_F`. |
| Unitarity | `matrixPrep` always produces a unitary. |
| Three-backend match | MS = MPS = SV for `matrixPrep`. |
| Row prep correctness | `buildRowQOp v |0..0⟩ = v / ‖v‖`. |

## Bugs caught so far

The current property suite has surfaced (and led to fixes for):

1. `permutationSwaps`'s `pos` array initialised as `p0` instead of `p0⁻¹`
   (only manifested on non-self-inverse permutations like 3-cycles).
2. `R axis θ` single-qubit case applied the bare Pauli, ignoring θ.
3. `R` multi-qubit sign error (`+i·s·P` instead of `−i·s·P`).
4. `SX` matrix `(p,p,m,p)` — non-unitary typo.
5. `dagger SX = Adjoint SX` infinite recursion.
6. `supportInterval` empty-support fallback out-of-bounds.
7. `MatrixPreparation.calculate{Real,Complex}Gate` sign convention.
8. `createRotationsDS` left-leaning DirectSum tree (dim-mismatch on
   non-power-of-2 chunk counts).
9. `op_support` for `Compose (Permute ks) b` returning ∅ when `b` had
   empty support.

## Backlog of invariants to add

Ordered roughly by expected bug-catch yield. Each item is a one-liner
property that can be encoded as either an HUnit case or a QuickCheck
property over `RandomQOp`.

### High-value, easy to write

**Simplifier preservation** — every rewrite rule in `HQP.QOp.Simplify`
should preserve `evalOp`. Untested today.

- `evalOp (cleanOnes op) ≈ evalOp op`
- `evalOp (cleanAdjoints op) ≈ evalOp op`
- `evalOp (liftComposes op) ≈ evalOp op`
- `evalOp (pushComposes op) ≈ evalOp op`
- `evalOp (doComposes op) ≈ evalOp op`
- `evalOp (cleanop op) ≈ evalOp op` (fixpoint composition)
- `evalOp (simplifyFixpoint n rules op) ≈ evalOp op` for any rule list.

**Adjoint distribution rules** — `dagger` is one of the most bug-prone
functions in the library; tested today only at the end-to-end level.

- `adj (Compose a b) = Compose (adj b) (adj a)` (note: reverses order)
- `adj (Tensor a b) = Tensor (adj a) (adj b)`
- `adj (DirectSum a b) = DirectSum (adj a) (adj b)`
- `adj (C a) = C (adj a)`
- `adj (R ax θ) = R ax (-θ)` (for Pauli `ax`)
- `adj (Permute π) = Permute (invertPerm π)`

**Algebraic structure** — operator-algebra laws. Catches whole classes
of backend evaluator bugs.

- `Compose` associativity: `(a ∘ b) ∘ c = a ∘ (b ∘ c)`
- `Tensor` associativity
- Identity laws: `a ∘ Id n = a`, `Id n ∘ a = a`, `a ⊗ One = a`,
  `One ⊗ a = a`
- Bifunctoriality: `(a ⊗ b) ∘ (c ⊗ d) = (a ∘ c) ⊗ (b ∘ d)` when arities
  align
- DirectSum-Compose distributivity: `(a ⊕ b) ∘ (c ⊕ d) = (a ∘ c) ⊕ (b ∘ d)`
- `C a = Id (n_a) ⊕ a` (cheap cross-backend check, since MS literally
  implements C this way)

### Medium-value, specific gate identities

**Pauli / Clifford algebra**

- `X² = Y² = Z² = H² = I`
- `(C X)² = I`
- `H ∘ X ∘ H = Z`, `H ∘ Z ∘ H = X`
- `(I ⊗ H) ∘ (C X) ∘ (I ⊗ H) = C Z`
- `R Z θ ∘ X = X ∘ R Z (-θ)`

**Permutation algebra** — likely to catch `applyPermute` /
`permutationSwaps` regressions.

- `Permute (compose π σ) = Permute π ∘ Permute σ`
- `Permute (identity perm) = Id n`
- `Permute π ∘ (X ⊗ Id (n-1)) ∘ Permute π⁻¹ = X` on qubit `π[0]`

**Rotation algebra** beyond what's there

- `R ax 4 = Id n_ax` (4π rotation)
- `R ax (a+b) = R ax a ∘ R ax b` (general angle additivity)
- For 1-qubit `ax`: `Tr(R ax θ) = 2 cos(πθ/2)`

**Phase**

- Matrix of `Phase a ∘ Phase b` = matrix of `Phase (a+b)`
- `Phase 0` = scalar 1, `Phase 2` = scalar 1

### Numerical / structural

**Norm and trace preservation**

- `inner (apply U ψ) (apply U ψ) = inner ψ ψ` for unitary U
- `Tr(U M U†) = Tr(M)`
- `|det U| = 1` for unitary U

**Under `Truncate` cfg in MPS**

- Bond dim never exceeds `maxBond`
- Reported `errorBound` in the profile is non-decreasing
- For Exact cfg with random ops, results are identical to MS within FP
  tolerance (already partially covered)

### Apply-path (state) properties

Useful once we want to test n > 4 qubits where full-matrix
materialisation is expensive.

- `apply (Tensor a b) (ket bs) = (apply a (ket b_a)) ⊗ (apply b (ket b_b))`
- `apply (Compose a b) ψ = apply a (apply b ψ)`
- Probabilities sum to 1 after a Born-rule measurement.

## Adding a property

1. Decide where it lives: cross-backend / general → `BackendValidation`;
   matrix-prep specific → `MatrixPrepTest`; algebraic identities →
   `AlgebraTest` (next file to add).
2. Write the property as `forAll genQOpAt …` or against `RandomQOp` /
   `RandomMatrix` / `RandomVec` from the generators above.
3. Match matrices via the `(~~)` operator (L1 diff < `tol`).
4. Register in the test exe's `defaultMain` group.
