{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE BangPatterns #-}

module HQP.QOp.MPSSemantics
  ( Site(..), MPS(..), StateT, WorkT, OpT(..)
  , Interval(..)
  , Trunc(..), EvalCfg(..), defaultCfg, ProfileCfg(..)
  , ket, ketW, toSparseMat, mpsToDenseVec, opToDenseMat
  , measureProjection, measure1, sampleAll
  , isDiagonalSite, isDiagonalMPS, siteVec
  , sampleAllDiagVecs, measureAllDiagVecs
  , sampleAllDiag, measureAllDiag
  , apply, evalStep, evalProg
  , evalOp, evalOpAtW
  , dagger
  , approxCfg, compressRange, moveCenterToPhys, bondDimensions, maxBondDimension
  ) where
-- | TODO: DONE Query bond dimension
---        DONE Track accumulated error
---        1. Implement "physical" permutation of MPS sites
---        2. R (1-qubit) θ -> single-site gate.
---        3. Clean up dirty bit (is it actually needed?)
---        4. Write documentation
---        5. Optimized multi-qubit rotation and multi-controlled gates
---        6. Drop WorkT, since we no longer have lazy evaluation?
---        7. Visualization tools for MPS

import Data.Complex (Complex(..), conjugate, magnitude, realPart, imagPart)
import qualified Data.Vector as V
import Data.Vector ((!), Vector, (//))
import Data.Bits (shiftL, testBit, (.|.),xor)
import Data.List (foldl')
import Data.Array (accumArray, elems)
import qualified Data.Set as S
import HQP.QOp.Syntax
import HQP.QOp.HelperFunctions
import HQP.PrettyPrint.PrettyMatrix
import HQP.PrettyPrint.PrettyOp
import HQP.QOp.MatrixSemantics (CMat, CVec)
import qualified HQP.QOp.MatrixSemantics as MS
import qualified Data.PQueue.Prio.Min as PriorityQ

import Numeric.LinearAlgebra (
  (><), (<#), (#>), atIndex, tr, rows, cols, size, diagBlock, diag, compactSVD, compactSVDTol, dot,
  takeRows,takeColumns,dropRows,dropColumns,cmap,maxElement, fromLists, toLists
  )
import qualified Numeric.LinearAlgebra as H
import Numeric.IEEE(epsilon)
import Numeric.LinearAlgebra.Devel (foldVector)
import Debug.Trace (trace)
import GHC.Stack (HasCallStack, callStack, prettyCallStack)


type StateT = MPS
type WorkT  = MPS
data OpT    = OpT { opQubits :: !Int, runOp :: MPS -> MPS }

apply :: OpT -> StateT -> StateT
apply (OpT _ f) x = f x

-- MPS representation
data Interval = Ival !Int !Int deriving (Show,Eq)
hull :: Interval -> Interval -> Interval
hull (Ival l r) (Ival l' r') = Ival (min l l') (max r r')

singletonIval :: Int -> Interval
singletonIval p = Ival p p

-- A : Dl x 2 x Dr as (A0,A1), each Dl x Dr.
data Site = Site { a0 :: !CMat, a1 :: !CMat } deriving (Show,Eq)
-- MPS in physical chain order; Permute updates only maps.
data MPS = MPS
  { scalar     :: !ComplexT
  , sites      :: !(Vector Site)
  , center_site :: !Int
  , log2phys   :: !(Vector Int)
  , phys2log   :: !(Vector Int)
  , dirty      :: !(Maybe Interval)
  , cfg        :: !EvalCfg

  } deriving (Show,Eq)

nSites :: MPS -> Int
nSites = V.length . sites

instance HasQubits MPS where n_qubits = nSites
instance HasQubits OpT where n_qubits = opQubits

class HasWork t where
  toWork   :: t -> WorkT
  fromWork :: WorkT -> t
instance HasWork StateT where toWork = id; fromWork = id

invertVec :: Vector Int -> Vector Int
invertVec v =
  let n = V.length v
  in  (V.replicate n 0) // [ (p,q) | (q,p) <- zip [0..] (V.toList v) ]

instance HasTensorProduct MPS where (⊗) = tensorMPS

{-| tensor product a ⊗ b for MPS a,b -}
tensorMPS :: MPS -> MPS -> MPS
tensorMPS a b =
  let (na,nb) = (nSites a, nSites b)
      (a',b') = (moveCenterToPhys (na-1) a, moveCenterToPhys 0 b)
      l2p = V.generate (na+nb) $ \q -> if q < na then log2phys a ! q else na + log2phys b ! (q-na)
  in MPS { scalar      = scalar a * scalar b,
           sites       = sites a' V.++ sites b',
           center_site = if na > 0 then center_site a else na + center_site b,
           log2phys = l2p,
           phys2log = invertVec l2p,
           dirty = Nothing,
           --dirty = Just (Ival (na-1) (na+1)),
           cfg = truncMax (cfg a) (cfg b) -- Worst accuracy guarantee determines overall accuracy
           }

{-| Combined truncation level: If we combine two MPS, use the least accurate truncation level -}
truncMax :: EvalCfg -> EvalCfg -> EvalCfg
truncMax c1 c2 = EvalCfg
  { trunc = case (trunc c1, trunc c2) of
      (Exact, t) -> t
      (t, Exact) -> t
      (Truncate m1 e1 p1, Truncate m2 e2 p2) -> Truncate (min m1 m2) (max e1 e2) (mergeProfile p1 p2)
  , tol = max (tol c1) (tol c2)
  }

mergeProfile :: Maybe ProfileCfg -> Maybe ProfileCfg -> Maybe ProfileCfg
mergeProfile Nothing p = p
mergeProfile p Nothing = p
mergeProfile (Just p1) (Just p2) = Just $ ProfileCfg
  { nSVDs      = nSVDs p1 + nSVDs p2
  , maxBondDim = max (maxBondDim p1) (maxBondDim p2)
  , errorBound = max (errorBound p1) (errorBound p2)
  }

withSameFrame :: String -> WorkT -> WorkT -> (WorkT -> WorkT -> a) -> a
withSameFrame ctx x y f
  | nSites x /= nSites y      = error (ctx ++ ": arity mismatch")
  | log2phys x == log2phys y  = f x y
  | otherwise                 = f (normalizeSiteOrder x) (normalizeSiteOrder y)

-- Relabel y so that its logical wires match the target frame.

scaleMPS :: ComplexT -> MPS -> MPS
scaleMPS c m = m { scalar = c * scalar m }

absorbScalarAt :: Int -> WorkT -> WorkT
absorbScalarAt p m
  | scalar m == (1:+0) = m
  | p < 0 || p >= nSites m = error "absorbScalarAt"
  | otherwise =
      let c = scalar m
          s = sites m ! p
          s' = Site (c .* (a0 s)) (c .* (a1 s))
      in m { scalar = 1:+0, sites = sites m // [(p,s')] }

markDirtyIval :: Interval -> MPS -> MPS
markDirtyIval i m = m { dirty = Just $ maybe i (hull i) (dirty m) }

clearDirty :: MPS -> MPS
clearDirty m = m { dirty = Nothing }

(.*.) :: CMat -> CMat -> CMat
(.*.) = (H.<>) -- HMatrix <> conflicts with Prelude.<>

innerMPS :: HasCallStack => WorkT -> WorkT -> ComplexT
innerMPS x0 y0 = -- trace ("innerMPS(" ++ showState x0 ++ ", " ++ showState y0 ++ ")") $
  withSameFrame "innerMPS" x0 y0 $ \x y ->
  let n  = nSites x
      (sx, sy) = (scalar x, scalar y)
      e0 = foldl' step ((1><1) [1:+0]) [n-1, n-2 .. 0] -- Full sweep for contraction
      step e p = -- e is contraction accumulator ("environment" in MPS lingo), p is current site
        let Site xa0 xa1 = sites x ! p
            Site ya0 ya1 = sites y ! p
            --dims xs = map (\m -> (rows m, cols m)) xs
            term xs ys = --trace(show $ dims [xs,tr xs,e,ys]) $
                        ys .*. e .*. tr xs -- Contract with accumulator. Why are dims reversed at this point?
        in term xa0 ya0 + term xa1 ya1 -- Contract over 0,1
      dotprod = conjugate sx * sy * (e0 `atIndex` (0,0))
  in
    --trace("innerMPS = " ++ show dotprod)
    dotprod

instance HilbertSpace WorkT where
  type Realnum WorkT = Double
  type Scalar  WorkT = ComplexT

  (.*) c ψ = scaleMPS c ψ
  (.+) ψ φ = addMPS ψ φ
  (.-) ψ φ = addMPS ψ (((-1):+0) .* φ)

  inner = innerMPS
  normalize ψ = let
     nrm = norm ψ
    in
      if nrm < tol (cfg ψ) then ψ else ((1/nrm):+0) .* ψ

data Trunc = Exact | Truncate { maxBond :: !Int, svd_r :: !Double, profile :: Maybe ProfileCfg }
        deriving (Show,Eq)

data EvalCfg    = EvalCfg { trunc :: !Trunc, tol :: !Double } deriving (Show,Eq)
data ProfileCfg = ProfileCfg {
  nSVDs      :: !Int,
  maxBondDim :: !Int,
  errorBound :: !Double
} deriving (Show,Eq)

defaultCfg :: EvalCfg
defaultCfg = EvalCfg { trunc = Exact, tol = 1e-12 }

approxCfg :: Bool -> Int -> Double -> EvalCfg
approxCfg doProfile maxbond svd_r = EvalCfg
  { trunc = Truncate { maxBond = maxbond, svd_r = svd_r,
                       profile = if doProfile then Just (ProfileCfg 0 0 0) else Nothing }
  , tol = 1e-12
  }

-- structural dagger -- move to Syntax?
dagger :: QOp -> QOp
dagger = \case
  Id n          -> Id n
  Phase q       -> Phase (-q)
  X             -> X
  Y             -> Y
  Z             -> Z
  H             -> H
  SX            -> Compose SX X       -- SX² = X, so SX⁻¹ = SX³ = SX·X

  R a t         -> R a (-t)
  C a           -> C (dagger a)
  Permute ks    -> Permute (invertPerm ks)
  Tensor a b    -> Tensor (dagger a) (dagger b)
  DirectSum a b -> DirectSum (dagger a) (dagger b)
  Compose a b   -> Compose (dagger b) (dagger a)
  Adjoint a     -> a

-- basis kets (logical MSB-first)
ketW :: [Int] -> WorkT
ketW bs =
  let n = length bs
      mk b = let v0 = if b==0 then 1:+0 else 0:+0
             in Site ((1><1) [v0]) ((1><1) [1-v0])
      ss  = V.fromList (map mk bs)
      idm = V.generate n id
  in MPS { scalar = 1:+0, sites = ss, center_site = 0, log2phys = idm, phys2log = idm, dirty = Nothing, cfg = defaultCfg }
ket :: [Int] -> StateT
ket = fromWork . ketW

-- wire permutation: logical relabeling only
applyPermute :: HasCallStack => Int -> [Int] -> WorkT -> WorkT
applyPermute base π m =
  let n    = length π
      l2p  = log2phys m
      sl   = V.slice base n l2p
      sl'  = V.fromList [ sl ! k | k <- π ]
      l2p' = l2p // [ (base+i, sl' ! i) | i <- [0..n-1] ]
      p2l' = invertVec l2p'
  in m { log2phys = l2p',
         phys2log = p2l' }

-- small matrix helpers
hcat, vcat :: CMat -> CMat -> CMat
hcat = (H.|||)
vcat = (H.===)

{-| Dense matrix representation of two adjacent sites
               [ A0 ]
     Θ2(A,B) = [ A1 ] [B0 B1] -}
theta2 :: Site -> Site -> CMat
theta2 a b = vcat (a0 a)  (a1 a) .*. hcat (a0 b) (a1 b)

diagMulLeft :: H.Vector Double -> CMat -> CMat
diagMulLeft v m = (H.complex . H.asColumn $ v) * m

diagMulRight :: CMat -> H.Vector Double -> CMat
diagMulRight m v = m * (H.complex . H.asRow $ v)

svd_compact :: EvalCfg -> CMat -> (CMat, H.Vector Double, CMat, Double)
svd_compact cfg m = case trunc cfg of
  Exact -> let
              (u,s,v) = compactSVD m
            in (u, s, v, 0)

  Truncate{maxBond,svd_r, profile} ->
    let
      (u,s,v) = svd profile m
      [u',v'] = H.takeColumns chi <$> [u,v] -- Truncate internal bond dimension to χ
      error_bound = case profile of
        Just _  -> H.sumElements (H.subVector chi (size s - chi) (s*s))
        Nothing -> 0

-- If we track errors, we need to calculate all nonzero singular values
      svd (Just _) = compactSVD ; svd Nothing = compactSVDTol svd_r
      svd_tol      = svd_r*g*epsilon*k where g = H.norm_Inf s
                                             k = fromIntegral (max (rows m) (cols m))
      chi          = min maxBond (firstBelow svd_tol s)
    in (u', s, v', error_bound)



updateProfile :: EvalCfg -> Int -- χ (new bond dim)
                -> Double         -- error ( ∑ σ_i^2 for discarded singular values σ)
                -> EvalCfg
updateProfile cfg chi error_bound = case trunc cfg of
  Truncate{profile=Just p}  ->
    let p' = p { nSVDs      = (nSVDs p) + 1,
                 maxBondDim = max (maxBondDim p) chi,
                 errorBound = errorBound p + error_bound
               }
    in cfg { trunc = (trunc cfg) { profile = Just p' } }
  _ -> cfg

moveRight :: Int -> WorkT -> WorkT
moveRight j m
  | j < 0 || j+1 >= nSites m = error "moveRight"
  | otherwise =
      let
          (a, b)        = (sites m ! j, sites m ! (j+1))
          (dl, dr)      = (rows (a0 a), cols (a0 b)) -- External bond dimensions
          (u,s,v,δ) = svd_compact (cfg m) (theta2 a b) -- compactSVD (exact or tol truncated)

          -- Update profiling statistics + track the error bound
          cfg' = updateProfile (cfg m) (size s) δ

          -- Absorb singular values into B'. Θ = U S V† = A' B'  =>  B' = S V†
          sv   = diagMulLeft s (tr v)
          a'   = uncurry Site (split2x1 dl u) -- Left Isometry A'
          b'   = uncurry Site (split1x2 dr sv)

      in m { sites = sites m // [(j,a'),(j+1,b')], center_site = j+1, cfg = cfg' }
-- factor to moveCenter
moveLeft :: Int -> WorkT -> WorkT
moveLeft j m
  | j <= 0 || j >= nSites m = error "moveLeft"
  | otherwise =
      let i = j-1
          (a,b)     = (sites m ! i, sites m ! j)
          (dl, dr)  = (rows (a0 a), cols (a0 b))
          (u,s,v,δ) = svd_compact (cfg m) (theta2 a b)

          -- Update profiling statistics + track the error bound
          cfg' = updateProfile (cfg m) (size s) δ

          -- Absorb singular values into A'. Θ = U S V† = A' B'  =>  A' = U S
          us  = diagMulRight u s
          a'  = uncurry Site (split2x1 dl us)
          b'  = uncurry Site (split1x2 dr (tr v))
      in m { sites = sites m // [(i,a'),(j,b')], center_site = j-1, cfg = cfg' }


-- TODO: This can be done cheaper by QR-decomposition and setting the dirty bit.
--       If inv(p) is quadratic, we can reduce to O(n^2) QR's and O(n) SVD's.
swapSites :: WorkT -> Int -> WorkT
swapSites psi j  = let
    n = nSites psi
  in
    if j < 0 || j+1 >= n then error "swapSites"
    else
      let (a,  b)  = (sites psi ! j, sites psi ! (j+1))
          (dl, dr) = (rows (a0 a), cols (a0 b))

          (a0b0,a0b1,
           a1b0,a1b1) = split2x2 dl dr $ theta2 a b

          theta' = H.fromBlocks [[a0b0,a1b0],
                                 [a0b1,a1b1]] -- Swap off-diagonal blocks

          (u,s,v,δ) = svd_compact (cfg psi) theta'

          a' = uncurry Site (split2x1 dl u)
          sv = diagMulLeft s (tr v)
          b' = uncurry Site (split1x2 dr sv)

          cfg' = updateProfile (cfg psi) (size s) δ
      in psi { sites = sites psi // [(j,a'),(j+1,b')], cfg = cfg' }

permutePhysicalSwaps :: WorkT -> [Int] -> WorkT
permutePhysicalSwaps = foldl' swapSites -- inv(pi) swaps (w/ SVD), so O(n^2 χ^3) worst case.

-- | Reorder the "physical" qubit sites to match the logical qubit order. This is needed before
--   adding two MPS together.
normalizeSiteOrder :: WorkT -> WorkT
normalizeSiteOrder psi =
  let swaps = permutationSwaps (phys2log psi)
      psi'  = permutePhysicalSwaps psi swaps
      ident = V.generate (nSites psi) id
  in psi' { log2phys = ident, phys2log = ident }

moveCenterToPhys :: Int -> WorkT -> WorkT
moveCenterToPhys p0 = go where
  go !m | center_site m < p0 = go (moveRight (center_site m) m)
        | center_site m > p0 = go (moveLeft  (center_site m) m)
        | otherwise         = m

compressRange :: Interval -> WorkT -> WorkT
compressRange (Ival l r) m
  | r <= l    = m
  | otherwise = foldl' (\acc j -> moveRight j acc) m [l..r-1]

compressIfDirty :: WorkT -> WorkT
compressIfDirty m = case dirty m of
  Nothing -> m
  Just i  -> clearDirty (compressRange i m)
-- local 1-qubit gate on a physical site: A'_s = sum_t U_{s,t} A_t
type Gate1 = (ComplexT,ComplexT,ComplexT,ComplexT) -- (u00,u01,u10,u11)

apply1Phys :: Gate1 -> Int -> WorkT -> WorkT
apply1Phys (u00,u01,u10,u11) p m =
  let s = sites m ! p
      a0' = u00 .* (a0 s) + u01 .* (a1 s)
      a1' = u10 .* (a0 s) + u11 .* (a1 s)
  in markDirtyIval (singletonIval p) $ m { sites = sites m // [(p, Site a0' a1')] }

apply1Logical :: Int -> Int -> Gate1 -> WorkT -> WorkT
apply1Logical base k u m = apply1Phys u (log2phys m ! (base+k)) m

pauliGate :: QOp -> Gate1
pauliGate = \case
  X -> (0,1,1,0)
  Y -> (0, 0:+(-1), 0:+1, 0)
  Z -> (1,0,0,-1)
  _ -> error "pauliGate"

cisPi :: Rational -> ComplexT
cisPi q = let t = pi * fromRational q in cos t :+ sin t
-- axis is Phase * Pauli-string in Tensor form (type-checked upstream)
-- returns (global phase phi, per-qubit operator list)
axisPaulis :: Int -> QOp -> (ComplexT, [QOp])
axisPaulis n = go (1:+0) where
  go !ϕ op = case op of
    Id m      -> (ϕ, replicate m (Id 1))
    Phase q   -> (ϕ * cisPi q, replicate n (Id 1))
    X         -> (ϕ, [X])
    Y         -> (ϕ, [Y])
    Z         -> (ϕ, [Z])

    Tensor a b ->
      let (ϕ1,p1) = axisPaulis (op_qubits a) a
          (ϕ2,p2) = axisPaulis (op_qubits b) b
      in (ϕ*ϕ1*ϕ2, p1++p2)

    -- TODO: Allow Compose


    Adjoint a -> let (ϕ1,ps) = axisPaulis n a in (conjugate ϕ * ϕ1, ps)

    _         -> error "axisPaulis: axis not Pauli string"

applyPauliString :: Int -> QOp -> WorkT -> (ComplexT, WorkT, Interval)
applyPauliString base axis m =
  let n = op_qubits axis
      (ϕ, ops) = axisPaulis n axis
      actIdx = op_support axis
      phys   = [ log2phys m ! (base+i) | i <- S.toList actIdx ]
      iSupp = case phys of
                [] -> singletonIval (log2phys m ! base)
                _  -> Ival (minimum phys) (maximum phys)
      m' = foldl' (\acc i -> case ops !! i of
                               Id _ -> acc
                               g    -> apply1Logical base i (pauliGate g) acc
                   ) m [0..n-1]
  in (ϕ, m', iSupp)

-- | Physical hull of positions where `phys2log[p] /= p`: exactly the sites
--   `normalizeSiteOrder` can mutate via `swapSites`, since the bubble sort
--   only sweeps between the lowest and highest misplaced positions. Reflects
--   ALL non-canonical state accumulated so far (e.g. log2phys inherited from
--   an outer `Permute`), so it can extend beyond any single op's
--   `op_support`. Returns `Nothing` when the state is already canonical.
normalizeTouchInterval :: WorkT -> Maybe Interval
normalizeTouchInterval m =
  let p2l = phys2log m
      n   = nSites m
      bad = [ p | p <- [0 .. n - 1], p2l ! p /= p ]
  in case bad of
       [] -> Nothing
       _  -> Just (Ival (minimum bad) (maximum bad))

-- local addition with branch-index only on a hull, followed by local compression
addLocal :: Interval -> WorkT -> WorkT -> WorkT
addLocal (Ival l0 r0) ψ0 φ0 = withSameFrame "addLocal" ψ0 φ0 $ \ψ φ ->
    -- The caller's `[l0, r0]` is derived from `op_support` — sites the current
    -- op touches. But when ψ0 and φ0 disagree on `log2phys`, `withSameFrame`
    -- invokes `normalizeSiteOrder`, whose `swapSites` (SVD-based) bubble sort
    -- sweeps the FULL range of phys2log misplacement — including misplacement
    -- inherited from ancestor ops, which can extend beyond any single op's
    -- support. After normalize, ψ and φ may disagree at sites the caller
    -- never named, so the merge's "outside [l, r] = identical" promise breaks.
    -- Fix: union [l0, r0] with each input's normalize-touch interval — the
    -- exact sweep range, no wider. Zero-cost when both inputs are canonical.
    let Ival l r =
          let combine Nothing iv = iv
              combine (Just (Ival a b)) (Ival a' b') = Ival (min a a') (max b b')
          in combine (normalizeTouchInterval ψ0)
                     (combine (normalizeTouchInterval φ0) (Ival l0 r0))
    in
    if      scalar ψ == 0 then φ0
    else if scalar φ == 0 then ψ0
    else let
        (ψ1,φ1) = (absorbScalarAt l ψ, absorbScalarAt l φ)
        (sψ, sφ) = (sites ψ1, sites φ1)

        -- Inside [l, r], combine ψ.site and φ.site by a position-dependent
        -- combinator: singleton support adds; left edge concatenates
        -- horizontally (entering Dl×2Dr); right edge concatenates vertically
        -- (exiting 2Dl×Dr); internal sites form a 2×2 block diagonal.
        -- Outside [l, r] the caller promised ψ.site p == φ.site p, so we
        -- take ψ.
        mk p
          | p < l || p > r = sψ ! p
          | otherwise =
              let a = sψ ! p; b = sφ ! p
                  combine = case (compare p l, compare p r) of
                    (EQ, EQ) -> (+)                      -- singleton support
                    (EQ, _ ) -> hcat                     -- left edge
                    (_,  EQ) -> vcat                     -- right edge
                    _        -> \x y -> diagBlock [x, y] -- internal
              in Site (combine (a0 a) (a0 b)) (combine (a1 a) (a1 b))

        out = ψ1 { scalar = 1:+0,
                   sites  = V.generate (nSites ψ1) mk }
  in compressRange (Ival l r) out

addMPS :: WorkT -> WorkT -> WorkT
addMPS ψ φ =
  let n = nSites ψ
      i = case (dirty ψ, dirty φ) of
            (Just a, Just b) -> hull a b
            (Just a, Nothing) -> a
            (Nothing, Just b) -> b
            (Nothing, Nothing) -> Ival 0 (max 0 (n-1))
  in addLocal i ψ φ

-- beam-search sparse extraction (fast for low-entanglement states)
toSparseMat :: (HasWork t) => Double -> Int -> t -> SparseMat
toSparseMat eps maxTerms t0 =
  let mps  = moveCenterToPhys 0 (toWork t0) -- ensures partial path amplitudes are strict bounds (yielding exact largest amplitudes)
      n    = nSites mps
      dim  = 2^(fromIntegral n :: Integer)
      eps2 = eps * eps

      step :: [(Integer, CVec)] -> Int -> [(Integer, CVec)]
      step states p =
        let Site a0 a1 = sites mps ! p
            q          = phys2log mps ! p
            bitpos     = (n - 1) - q
            branches   = [(0, a0), (1, a1)] :: [(Integer,CMat)]

            extend1
              :: (Integer, CMat)
              -> (Integer, CVec)
              -> PriorityQ.MinPQueue Double (Integer, CVec)
              -> PriorityQ.MinPQueue Double (Integer, CVec)
            extend1 (s,a) (!idx,!v) beam0 =
              let !v' = v <# a
                  !w2 = realPart $ dot v' v'
              in  if w2 <= eps2
                    then beam0 -- Continuing on this path cannot yield amplitude larger than eps2
                    else let !idx' = idx .|. (s `shiftL` bitpos)   -- Feasible candidate:
                         in  pushTopK maxTerms w2 (idx', v') beam0 -- push to top-k priority queue

            extendState
              :: PriorityQ.MinPQueue Double (Integer, CVec)
              -> (Integer, CVec)
              -> PriorityQ.MinPQueue Double (Integer, CVec)
            extendState beam0 st = foldl' (\b br -> extend1 br st b) beam0 branches

            beam = foldl' extendState PriorityQ.empty states
        in  map snd (PriorityQ.toDescList beam)

      finals = -- Branch and bound on path through sites from left to right
        foldl' step [(0 :: Integer, H.fromList [scalar mps])] [0 .. n-1]

      nz = [ ((i,0), v `atIndex` 0) | (i,v) <- finals ]
  in SparseMat ((dim,1), nz)


instance Convertible WorkT SparseMat where
  to   mps = toSparseMat (tol $ cfg mps) 100 mps
  from = \sm ->
    let (SparseMat ((m,_), nonzeros)) = sm
        n = integerlog2 m
        kets = [ a .* ketW (toBits' n k) | ((k,_),a) <- nonzeros ] :: [WorkT]
    in case kets of
         [] -> error "fromSparseMat: empty"
         (x:xs) -> foldl' (.+) x xs


instance Convertible StateT CMat where
  to   = mpsToDenseVec
  from psimat =
    let psi_sparse = MS.sparseMat psimat
    in from psi_sparse :: StateT

-- | Contract an MPS into a dense (2^n × 1) column vector, indexed by logical bits MSB-first
--   (matching MatrixSemantics' ket convention). Faithfully reflects the MPS's stored
--   amplitudes — does not modify the MPS or add any truncation beyond what is already in it.
--   O(2^n · χ²); intended for small n.
mpsToDenseVec :: MPS -> CMat
mpsToDenseVec mps =
  let n = nSites mps
      step !t p =
        let Site x0 x1 = sites mps ! p
            t0 = t .*. H.tr' x0
            t1 = t .*. H.tr' x1
        in t0 H.=== t1
      vPhys = H.scale (scalar mps) $
                foldl' step ((1><1) [1:+0]) [n-1, n-2 .. 0]
      l2p   = log2phys mps
      logToPhys !iLog = foldl' (.|.) 0
        [ if testBit iLog (n-1-q) then 1 `shiftL` (n-1-(l2p ! q)) else 0
        | q <- [0..n-1] ]
      dim   = pow2 n
  in H.asColumn $ H.fromList
       [ vPhys `atIndex` (logToPhys iLog, 0) | iLog <- [0 .. dim-1] ]

-- | Materialize an `OpT` as a dense matrix by applying it to each canonical basis ket
--   prepared with the given EvalCfg, then assembling the resulting column vectors. The cfg
--   propagates into the lambda's internal MPS operations, so the matrix reflects the
--   operator's actual truncation behavior under that cfg. Intended for small n.
opToDenseMat :: EvalCfg -> OpT -> CMat
opToDenseMat cfg' (OpT n f) =
  let mkKet bits = (ketW bits) { cfg = cfg' }
      cols       = [ mpsToDenseVec (f (mkKet (toBits' n j))) | j <- [0 .. pow2 n - 1] ]
  in H.fromBlocks [cols]

instance Convertible OpT CMat where
  to     = opToDenseMat defaultCfg
  from _ = error "Convertible OpT CMat: from is not implemented; build OpT via evalOp"

-- measurement
frob2 :: CMat -> Double
frob2 x =
  let v = H.flatten x
  in realPart $ dot v v -- Hmatrix conjugates left argument to dot

-- | TODO: Clean up compressIfDirty mess.
projectCtrl :: Int -> Bool -> WorkT -> WorkT
projectCtrl p one m =
  let m0 = moveCenterToPhys p (compressIfDirty m)
  in projectCenter p one  m0

-- | Project the *center* site to |b> (b=False => |0>, True => |1>).
--   Produces an *unnormalized* post-measurement state.
projectCenter :: Int -> Bool -> MPS -> MPS
projectCenter p b st =
  let Site x0 x1 = sites st ! p
      (y0,y1)    = if b then (zeros_like x0, x1) else (x0, zeros_like x1)
  in st { sites = sites st // [(p, Site y0 y1)]
        , dirty = Just (singletonIval p)
        }

measureProjection :: HasCallStack => Int -> Int -> Int -> OpT
measureProjection arity k out = OpT arity $ \t0 ->
  let st0 = compressIfDirty t0
      p   = log2phys st0 ! k
      st2 = projectCenter p (out==1) st0
      st3 = clearDirty (compressRange (singletonIval p) st2)
  in st3



zeros_like :: CMat -> CMat
zeros_like x = H.konst (0:+0) (rows x, cols x)

-- | TODO: Simplify
-- | TODO: Clean up compressIfDirty mess.
measure1 :: HasCallStack => (StateT, Outcomes, RNG) -> Int -> (StateT, Outcomes, RNG)
measure1 (st, outs, u:us) k = --trace("measure1 on qubit " ++ show k ++ " of state " ++ showState st) $
  let
      st0 = compressIfDirty st
      p   = log2phys st0 ! k
      st1 = moveCenterToPhys p st0
      Site x0 x1 = sites st1 ! p
      s2  = let a = magnitude (scalar st1) in a*a
      p0  = s2 * frob2 x0
      p1  = s2 * frob2 x1
      tot = p0 + p1
  in if tot < (tol $ cfg st) then error "measure1: prob~0" else
     let b    = (u*tot >= p0)
         pb   = if b then p1 else p0
         inv  = (1 / sqrt pb) :+ 0
         y0   = if b then zeros_like x0 else inv .* x0
         y1   = if b then inv .* x1 else zeros_like x1
         st2  = st1 { sites = sites st1 // [(p, Site y0 y1)], dirty = Just (singletonIval p), cfg = cfg st1 }
         st3  = clearDirty (compressRange (singletonIval p) st2)
     in (st3, b:outs, us)
measure1 (_,_,[]) _ = error "measure1: empty RNG"

-- | Sample a computational-basis outcome for *every* qubit without producing a
--   post-measurement state. Uses the Ferris--Vidal perfect-sampling recurrence:
--   given a mixed-canonical MPS at center @c@, sweep right from @c@ to @n-1@
--   maintaining a running matrix product @B@ as the left boundary, then sweep
--   left from @c-1@ to @0@ where the boundary has collapsed to a vector. No
--   SVDs; matmul for the right sweep, matvec for the left sweep.
--   Outcomes are returned head-most-recent, matching @Measure [0..n-1]@:
--   @outs = [b_0, b_1, ..., b_{n-1}]@ with @b_0@ at the head (last bit
--   processed by the @reverse ks@ fold).
sampleAll :: HasCallStack => StateT -> RNG -> (Outcomes, RNG)
sampleAll = sampleAllW . toWork

sampleAllW :: HasCallStack => WorkT -> RNG -> (Outcomes, RNG)
sampleAllW psi0 rng0
  | n == 0    = ([], rng0)
  | otherwise =
      let chiC  = rows (a0 (sites psi ! c))
          bInit = H.complex (H.ident chiC :: H.Matrix Double) :: CMat
          (bEnd, rng1, rightPairs) =
            foldl' rightStep (bInit, rng0, []) [c .. n-1]
          vInit = H.flatten bEnd
          (_v,   rng2, leftPairs)  =
            foldl' leftStep  (vInit, rng1, []) [c-1, c-2 .. 0]
          bits = V.replicate n False // (rightPairs ++ leftPairs)
          outs = [ bits ! k | k <- [0 .. n-1] ]
      in (outs, rng2)
  where
    psi   = compressIfDirty psi0
    n     = nSites psi
    c     = center_site psi
    tolP  = tol (cfg psi)

    rightStep :: (CMat, RNG, [(Int,Bool)]) -> Int
              -> (CMat, RNG, [(Int,Bool)])
    rightStep (_, [], _)     _ = error "sampleAll: empty RNG"
    rightStep (b, r:rs, acc) p =
      let Site x0 x1 = sites psi ! p
          m0  = b .*. x0
          m1  = b .*. x1
          w0  = frob2 m0
          w1  = frob2 m1
          tot = w0 + w1
      in if tot < tolP
           then error "sampleAll: prob ~ 0"
           else let bit = r * tot >= w0
                    b'  = if bit then m1 else m0
                    q   = phys2log psi ! p
                in (b', rs, (q, bit) : acc)

    leftStep :: (CVec, RNG, [(Int,Bool)]) -> Int
             -> (CVec, RNG, [(Int,Bool)])
    leftStep (_, [], _)     _ = error "sampleAll: empty RNG"
    leftStep (v, r:rs, acc) p =
      let Site x0 x1 = sites psi ! p
          u0  = x0 #> v
          u1  = x1 #> v
          w0  = realPart (dot u0 u0)   -- dot conjugates left arg => ‖u‖²
          w1  = realPart (dot u1 u1)
          tot = w0 + w1
      in if tot < tolP
           then error "sampleAll: prob ~ 0"
           else let bit = r * tot >= w0
                    v'  = if bit then u1 else u0
                    q   = phys2log psi ! p
                in (v', rs, (q, bit) : acc)

------------------------------------------------------------------------
-- Diagonal-MPS specialization
------------------------------------------------------------------------
-- A "diagonal-core" MPS represents a sum of χ separable terms:
--     |ψ⟩ = Σ_α ⊗_p (f_p(0,α) |0⟩ + f_p(1,α) |1⟩)
-- In the standard MPS chain, this means boundary sites are 1×χ / χ×1
-- (row/column vector pairs) and internal sites are χ×χ diagonal in the
-- bond index. Sampling collapses to O(n·χ²) — see perfect-sampling.tex
-- section "Diagonal-MPS specialization".
--
-- Two entry-point layers:
--   * Raw (sampleAllDiagVecs / measureAllDiagVecs) — takes per-site
--     (f_p(0,·), f_p(1,·)) vector pairs directly. No χ² storage required
--     by the caller; only O(χ²) working memory for the R-chain.
--   * MPS wrappers (sampleAllDiag / measureAllDiag) — assert
--     isDiagonalMPS, extract via siteVec, run the raw entry, remap
--     outcomes to logical order via phys2log.

-- | Per-site vector pair (f_p(0,·), f_p(1,·)) extracted from a Site,
--   handling the three shape cases uniformly.
siteVec :: Site -> (CVec, CVec)
siteVec (Site x0 x1)
  | rows x0 == 1 = (H.flatten x0, H.flatten x1)
  | cols x0 == 1 = (H.flatten x0, H.flatten x1)
  | otherwise    = (H.takeDiag x0, H.takeDiag x1)

-- | A site is diagonal-or-boundary if it's 1×χ, χ×1, or square with all
--   off-diagonals below @eps@ in magnitude.
isDiagonalSite :: Double -> Site -> Bool
isDiagonalSite eps (Site x0 x1) =
  let isDiag m =
        let r = rows m; c = cols m
            offDiagSmall =
              and [ magnitude (m `atIndex` (i,j)) < eps
                  | i <- [0 .. r-1], j <- [0 .. c-1], i /= j ]
        in (r == 1) || (c == 1) || (r == c && offDiagSmall)
  in isDiag x0 && isDiag x1

-- | Whole-MPS predicate. Requires that boundary sites have the correct
--   shape and all internal sites are diagonal. Uses @tol (cfg psi)@.
isDiagonalMPS :: MPS -> Bool
isDiagonalMPS psi =
  let n      = nSites psi
      sV     = sites psi
      eps    = tol (cfg psi)
      sLeft  = sV ! 0
      sRight = sV ! (n-1)
  in case n of
       0 -> True
       1 -> True   -- 1×1
       _ -> rows (a0 sLeft) == 1
         && cols (a0 sRight) == 1
         && all (\p -> isDiagonalSite eps (sV ! p)) [1 .. n-2]

-- | G_p[α,β] = Σ_s f_p(s,α)* · f_p(s,β)  — Hermitian, rank ≤ 2.
--   (Note conjugation on the *first* index: this is what falls out of
--   |⟨s|ψ⟩|² when integrating out one site.)
gMatrix :: (CVec, CVec) -> CMat
gMatrix (v0, v1) =
  let outerC u = H.outer (cmap conjugate u) u
  in outerC v0 + outerC v1

-- | Right-chain @R[p][α,β] = ∏_{q > p} G_q[α,β]@, computed as elementwise
--   products. @R[n-1]@ is the all-ones χ×χ matrix.
buildRChain :: Int -> V.Vector (CVec, CVec) -> V.Vector CMat
buildRChain chi vecs =
  let n     = V.length vecs
      jOnes = H.konst (1:+0) (chi, chi)
      gs    = V.map gMatrix vecs    -- G_0 .. G_{n-1}
      -- R[p] = G[p+1] * R[p+1] elementwise (Matrix Num instance is elementwise).
      go p acc | p < 0     = acc
               | otherwise = let r' = if p == n-1 then jOnes
                                                  else (gs ! (p+1)) * (acc ! (p+1))
                             in go (p-1) (acc // [(p, r')])
  in go (n-1) (V.replicate n jOnes)

-- | Sample all qubits from a diagonal-MPS given only per-site vector pairs.
--   Outcomes are returned head-most-recent, matching @Measure [0..n-1]@.
--   Cost: O(n · χ²); working memory: O(n · χ²) for the R-chain.
sampleAllDiagVecs :: HasCallStack
                  => V.Vector (CVec, CVec)
                  -> RNG
                  -> (Outcomes, RNG)
sampleAllDiagVecs vecs rng0
  | n == 0    = ([], rng0)
  | otherwise =
      let chi    = H.size (fst (vecs ! 0))
          rChain = buildRChain chi vecs
          lInit  = H.konst (1:+0) chi
          (_, rng', bitsRev) =
            foldl' (step rChain) (lInit, rng0, []) [0 .. n-1]
      in (reverse bitsRev, rng')
  where
    n = V.length vecs
    -- one Born-rule decision per site; uses chi from L's current size.
    step :: V.Vector CMat
         -> (CVec, RNG, [Bool])
         -> Int
         -> (CVec, RNG, [Bool])
    step _ (_, [],   _)   _ = error "sampleAllDiagVecs: empty RNG"
    step rChain (l, r:rs, acc) p =
      let (f0, f1) = vecs ! p
          h0  = l * f0           -- elementwise (Vector Num is elementwise)
          h1  = l * f1
          r_p = rChain ! p
          w0  = realPart (dot h0 (r_p #> h0))
          w1  = realPart (dot h1 (r_p #> h1))
          tot = w0 + w1
      in if tot < 1e-300
           then error "sampleAllDiagVecs: prob ~ 0"
           else let bit = r * tot >= w0
                    l'  = if bit then h1 else h0
                in (l', rs, bit : acc)

-- | Projective measure-all on a diagonal-MPS given vector pairs. Returns:
--     * post-measurement vectors (unchosen branch zeroed per site),
--     * renormalization factor 1/|⟨s|ψ⟩| (multiply into MPS scalar),
--     * outcomes (head-most-recent),
--     * remaining RNG.
--   The post-measurement vectors remain in the diagonal class.
measureAllDiagVecs :: HasCallStack
                   => V.Vector (CVec, CVec)
                   -> RNG
                   -> ( V.Vector (CVec, CVec)
                      , Complex Double
                      , Outcomes
                      , RNG )
measureAllDiagVecs vecs rng0
  | n == 0    = (vecs, 1:+0, [], rng0)
  | otherwise =
      let chi    = H.size (fst (vecs ! 0))
          rChain = buildRChain chi vecs
          zerof  = H.konst (0:+0) chi
          lInit  = H.konst (1:+0) chi
          (lFinal, rng', bitsRev, vecsOut) =
            foldl' (step rChain zerof) (lInit, rng0, [], vecs) [0 .. n-1]
          amp     = H.sumElements lFinal      -- ⟨s|ψ_internal⟩
          renorm  = (1 / magnitude amp) :+ 0
      in (vecsOut, renorm, reverse bitsRev, rng')
  where
    n = V.length vecs
    step :: V.Vector CMat
         -> CVec
         -> (CVec, RNG, [Bool], V.Vector (CVec, CVec))
         -> Int
         -> (CVec, RNG, [Bool], V.Vector (CVec, CVec))
    step _ _ (_, [], _, _) _ = error "measureAllDiagVecs: empty RNG"
    step rChain zerof (l, r:rs, acc, vs) p =
      let (f0, f1) = vs ! p
          h0  = l * f0
          h1  = l * f1
          r_p = rChain ! p
          w0  = realPart (dot h0 (r_p #> h0))
          w1  = realPart (dot h1 (r_p #> h1))
          tot = w0 + w1
      in if tot < 1e-300
           then error "measureAllDiagVecs: prob ~ 0"
           else let bit = r * tot >= w0
                    l'  = if bit then h1 else h0
                    vp' = if bit then (zerof, f1) else (f0, zerof)
                in (l', rs, bit : acc, vs // [(p, vp')])

-- | Sample-only all-qubit measurement specialized to diagonal MPS.
--   Asserts diagonality on entry. Outcomes returned in logical order.
sampleAllDiag :: HasCallStack => StateT -> RNG -> (Outcomes, RNG)
sampleAllDiag psi rng
  | not (isDiagonalMPS psi) = error "sampleAllDiag: state is not diagonal"
  | otherwise =
      let vecs           = V.map siteVec (sites psi)
          (physOuts, rng') = sampleAllDiagVecs vecs rng
          n              = nSites psi
          bits = V.replicate n False //
                   [ (phys2log psi ! p, b) | (p, b) <- zip [0..] physOuts ]
          outs = [ bits ! k | k <- [0 .. n-1] ]
      in (outs, rng')

-- | Projective all-qubit measurement on a diagonal MPS. Returns the
--   post-measurement state (still diagonal) with global scalar updated to
--   keep ψ normalised.
measureAllDiag :: HasCallStack => StateT -> RNG -> (StateT, Outcomes, RNG)
measureAllDiag psi rng
  | not (isDiagonalMPS psi) = error "measureAllDiag: state is not diagonal"
  | otherwise =
      let vecs                       = V.map siteVec (sites psi)
          (vecsOut, renorm, physOuts, rng') = measureAllDiagVecs vecs rng
          n        = nSites psi
          newSites = V.zipWith updateSiteFromVec (sites psi) vecsOut
          -- Post-projection: state magnitude is |old_scalar · ⟨s|ψ_internal⟩|.
          -- We want unit norm. The renorm factor from the raw entry is
          -- 1/|⟨s|ψ_internal⟩|, so setting scalar = renorm gives a state of
          -- magnitude 1 regardless of the original scalar (phase is
          -- unobservable post-collapse).
          psi'     = psi { sites = newSites
                         , scalar = renorm
                         }
          bits = V.replicate n False //
                   [ (phys2log psi ! p, b) | (p, b) <- zip [0..] physOuts ]
          outs = [ bits ! k | k <- [0 .. n-1] ]
      in (psi', outs, rng')

-- | Build a Site from a (CVec, CVec) given an original Site (for shape).
updateSiteFromVec :: Site -> (CVec, CVec) -> Site
updateSiteFromVec (Site origA0 _) (v0, v1)
  | rows origA0 == 1 = Site (H.asRow v0)    (H.asRow v1)
  | cols origA0 == 1 = Site (H.asColumn v0) (H.asColumn v1)
  | otherwise        = Site (H.diag v0)     (H.diag v1)

-- support interval (physical hull) from op_support
supportInterval :: WorkT -> Int -> QOp -> Interval
supportInterval st base op =
  let supL = S.toList (op_support op)
      ps   = [ log2phys st ! (base+q) | q <- supL ]
  in case ps of
       [] -> singletonIval (log2phys st ! base)
       _  -> Ival (minimum ps) (maximum ps)

-- evaluator
evalOp :: QOp -> OpT
evalOp op = OpT (op_qubits op) (evalOpAtW 0 op)


evalOpAtW :: HasCallStack => Int -> QOp -> WorkT -> WorkT
evalOpAtW base op st = case op of
  Id _      -> st
  Phase q   -> (cisPi q) .* st
  Permute π -> applyPermute base π st
  X  -> apply1Logical base 0 (pauliGate X) st
  Y  -> apply1Logical base 0 (pauliGate Y) st
  Z  -> apply1Logical base 0 (pauliGate Z) st
  H  -> let s = (1/sqrt 2):+0 in apply1Logical base 0 (s,s,s,-s) st
  SX ->
    let p = 0.5:+0.5
        m = 0.5:+(-0.5)
    in apply1Logical base 0 (p, m, m, p) st

  Tensor a b ->
    let st1 = evalOpAtW base a st
    in evalOpAtW (base + op_qubits a) b st1

  Compose a b ->
     evalOpAtW base a (evalOpAtW base b st)

  Adjoint a -> evalOpAtW base (dagger a) st

  R axis θ
    | θ == 0 -> st
    | otherwise ->
        let t = pi * fromRational θ / 2
            c = cos t :+ 0
            s = sin t :+ 0
            iC = 0 :+ 1
            (phi, pst, iSupp) = applyPauliString base axis st
            ψ1 = c .* st
            ψ2 = ((-iC)*s*phi) .* pst   -- exp(-iπθ/2 · P) = c·I − i·s·P
        in addLocal iSupp ψ1 ψ2

  C a ->
    let ctrlP  = log2phys st ! base
        ctrlIv = singletonIval ctrlP
        iA     = if S.null (op_support a)
                   then ctrlIv                          -- empty support: stay at control bit
                   else supportInterval st (base+1) a
        iHull  = hull ctrlIv iA
        psi0   = projectCtrl ctrlP False st
        psi1   = evalOpAtW  (base+1) a (projectCtrl ctrlP True st)
    in addLocal iHull psi0 psi1

  DirectSum a b ->
    let ctrlP  = log2phys st ! base
        ctrlIv = singletonIval ctrlP
        iA     = if S.null (op_support a) then ctrlIv else supportInterval st (base+1) a
        iB     = if S.null (op_support b) then ctrlIv else supportInterval st (base+1) b
        iHull  = hull ctrlIv (hull iA iB)
        psi0   = evalOpAtW (base+1) a (projectCtrl ctrlP False st)
        psi1   = evalOpAtW (base+1) b (projectCtrl ctrlP True  st)
    in addLocal iHull psi0 psi1

-- Steps / programs -- MOVE TO COMMON MODULE.
evalStep :: HasCallStack => (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalStep (st, outs, rng) step = -- trace("step "++showStep step ++ " on " ++ showState st) $
  case step of
  Unitary op -> (apply (evalOp op) st, outs, rng)
  Measure ks -> foldl' measure1 (st, outs, rng) (reverse ks)
  Initialize ks vs ->
      let
        n = n_qubits st
        (st', os, rng') = evalStep (st, [], rng) (Measure ks)
        -- List of outcomes xor values for each initialized qubit
        corrections     = zipWith xor os vs
        -- Now we build the full list, including unaffected qubits
        corrFull        = accumArray xor False (0,n-1) (zip ks corrections)
        corrOp          = foldl (⊗) One [ if c then X else I | c <- elems corrFull ]
      in (evalStep) (st', outs, rng') (Unitary corrOp)

evalProg :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
evalProg prog st rng = foldl' evalStep (st, [], rng) prog

-- Helper functions for MPS
bondDimensions :: MPS -> V.Vector Int
bondDimensions mps = bondDim <$> (sites mps)
  where
    bondDim (Site a0 _) = rows a0


maxBondDimension :: MPS -> Int
maxBondDimension mps = maximum . V.toList . bondDimensions $ mps
