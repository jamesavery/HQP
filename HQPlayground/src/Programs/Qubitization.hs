module Programs.Qubitization where
import HQP
import Polynomial.Roots (roots) -- From dsp package
import Polynomial.Basic (polyadd, polysub, polymult,polyderiv, polyeval)
import Data.Complex (Complex((:+)), realPart, imagPart,magnitude,polar)
import qualified Data.Vector as V
import Data.List (sortOn,sortBy, groupBy, partition)
import Data.Ord (comparing)
--import Debug.Trace

--------------------------------------------------------------------------------
------------- Quantum signal processing-----------------------------------
--------------------------------------------------------------------------------

-- qspVectorPhi calculates the quantum signal processing vector  
-- defined in Gilyen et al. Theorem 3.
qspVectorPhi :: [ComplexT] -> [ComplexT] -> [Double]
qspVectorPhi p q
    | polyDegree p == 0 =
        let 
            (_, phi) = polar (head p)
        in  
            [phi]           

    | otherwise =
        let
            -- Calculate phi value
            leadingOrderTermP = case getLeadingOrderTerm p of
                                Just (c, _deg) -> c
                                Nothing        -> error "getLeadingOrderTerm: Polynomial is empty or zero."

            leadingOrderTermQ = case getLeadingOrderTerm q of
                                Just (c, _deg) -> c
                                Nothing        -> error "getLeadingOrderTerm: Polynomial is empty or zero."

            lTermFrac = leadingOrderTermP / leadingOrderTermQ
            (_, theta) = polar lTermFrac
            phi = theta / 2

            xPoly =  [0,1]
            x2Poly = [1, 0, -1]
            pNextTerm1 = polymult xPoly p
            pNextTerm2 =  map (* lTermFrac) (polymult x2Poly q)
            pNext = map (* exp (-(0.0 :+ phi))) ( polyadd pNextTerm1 pNextTerm2 )

            qNextTerm1 = map (* lTermFrac) (polymult xPoly q)
            qNext = map (* exp (-(0.0 :+ phi))) ( polysub qNextTerm1 p)
        in
            (qspVectorPhi pNext qNext) ++ [phi]


--------------------------------------------------------------------------------
------------- From real to complex polynomial-----------------------------------
--------------------------------------------------------------------------------

-- Given a real polynomial satisfying the conditions:
-- ... KOMMER SENERE ... and
-- P(x)^2 <= 1 for x in [-1,1]
-- qspPolys calculates the pair of complex polynomials satisfying the conditions in Gilyen et al. (2018) Theorem 5.
-- The polynomials are used in algorithm (qspVectorPhi) that finds the quantum signal processing vector.
qspPolys :: [Double] -> ([ComplexT], [ComplexT])
qspPolys realEvenPosPoly 
    | not (checkPolyBounds realEvenPosPoly) = 
        error "createRealPolyPair: Condition P(x)^2 <= 1 not met in [-1,1] (Gilyen et al. (2018) Thm.5)"

    | otherwise = 
        let 
            sqPoly = polymult realEvenPosPoly realEvenPosPoly
            realPolyPair = createRealPolyPair (polysub [1]  sqPoly)
            
            realPoly1 = poly1 realPolyPair
            realPoly2 = poly2 realPolyPair

            -- Type convertions
            complexPoly1 =  map (:+ 0.0) realPoly1
            complexPoly2 =  map (:+ 0.0) realPoly2
            complexEvenPosPoly = map (:+ 0.0) realEvenPosPoly

            pPoly = polyadd complexEvenPosPoly (map (* (0.0 :+ 1.0)) complexPoly1)
            qPoly = map (* (0.0 :+ 1.0)) complexPoly2
        in
            (pPoly,qPoly)

-- Given a real even poly A(x) with deg(A) <= 2k and A(x)>= 0 on [-1,1] 
-- return real polys B, C as in Gilyen et al. Lemma 6.
-- The method is a step in the construction of a suitable complex polynomial
-- used in the qubitization algorithm
createRealPolyPair :: [Double] -> PolyPair
createRealPolyPair [] = error "createRealPolyPair: Coefficient list is empty."
createRealPolyPair coeffs = 
    let 
        epsilon = 1e-9
        maxIter = 1000

        -- 1. Base scale (needed regardless of roots)
        leadingCoeff = case getLeadingOrderTerm coeffs of
            Just (c, _)  -> c
            _            -> error "createRealPolyPair: Polynomium is zero."
        initialScale = PolyPair ([sqrt (abs leadingCoeff)],[0])

        -- 2. Find all roots
        allRoots = roots epsilon maxIter (map (:+ 0) coeffs)

    in case allRoots of
        -- EARLY EXIT: If no roots exist, just return the scale.
        [] -> initialScale
        
        -- PROCESS ROOTS: Only executes if allRoots has content.
        rs -> 
            let 
                -- Partition into [0, 1) real roots and others
                (subOne, others) = partition isSubOne rs
                    where 
                        isSubOne (r :+ i) = abs i < epsilon && r >= - epsilon && r < (1 - epsilon)

                -- Group 1: Sub-one roots (Cluster, check even multiplicity, process)
                subOnePairs = case clusterRoots 0.01 subOne of
                    [] -> []
                    cs -> if any (odd . snd) cs
                          then error $ "Multiplicity check failed near: " ++ show (map fst (filter (odd . snd) cs))
                          else processAllsubOneRoots cs

                -- Group 2: Remaining relevant roots (Filter and map)
                remainingPairs = map rootToPolyPair $ 
                    filter (\(r :+ i) -> (r >= 0 && i > 0) || (r >= 1 && i == 0)) others

            in foldr (*) initialScale (subOnePairs ++ remainingPairs)



-- Used to avoid multiplicity issues in the numerical root finding algorithm
clusterRoots :: Double -> [ComplexT] -> [(ComplexT, Int)]
clusterRoots tol = foldr integrate [] 
  where
    integrate r [] = [(r, 1)]
    integrate r ((center, count):rest)
        | magnitude (center - r) < tol = 
            let n = fromIntegral count
                newCenter = (center * n + r) / (n + 1)
            in (newCenter, count + 1) : rest
        | otherwise = (center, count) : integrate r rest


-- Handling of even multiplicity roots < 1
processAllsubOneRoots :: [(ComplexT, Int)] -> [PolyPair]
processAllsubOneRoots clusters =
        concatMap (\(root, mult) -> replicate (mult `div` 2) (rootToPolyPair root)) clusters

rootToPolyPair :: ComplexT -> PolyPair
rootToPolyPair z
    -- Case: Re z > 0 and Im z > 0
    | re >= epsilon && im >= epsilon = 
        let c = re^2 + im^2 + sqrt( (2 * (re^2 + 1) * im^2) + (re^2 - 1)^2 + im^4 )
            rePol = [-(re^2 + im^2), 0, c]
            imPol = [0, sqrt(c^2 - 1)]
        in 
            PolyPair (rePol, imPol)

    -- Case: z = 0
    | abs re < epsilon && abs im < epsilon = 
        let
            rePol = [0,1]
            imPol = [0]
        in
            PolyPair (rePol, imPol)

    -- Case: z in (0, 1) (Real axis between 0 and 1)
    | abs im < epsilon && re >= epsilon && re <= 1 - epsilon = 
        let
            rePol = [-re^2,0,1]
            imPol = [0]
        in
            PolyPair (rePol, imPol)

    -- Case: z in [1, infinity) (Real axis >= 1)
    | abs im < epsilon && re > 1 - epsilon = 
        let
            rePol = [0,sqrt(re^2 - 1)]
            imPol = [re]
        in
            PolyPair (rePol, imPol)

    -- Case: Re z = 0 and Im z > 0 (Purely imaginary)
    | abs re < epsilon && im >= epsilon = 
        let
            rePol = [0,sqrt(im^2 + 1)]
            imPol = [im]
        in
            PolyPair (rePol, imPol)

    -- Fallback/Default
    | otherwise = error "Unspecified complex region"
  where
    re = realPart z
    im = imagPart z
    epsilon = 1e-9

-- Polypair is a pair of real polyomials with a distint multipication
-- used in Gilyen et al. Lemma 6.
newtype PolyPair = PolyPair ([Double], [Double])
    deriving (Show, Eq)

instance Num PolyPair where
    -- (p, q) + (s, t) = (p+s, q+t)
    (PolyPair (p, q)) + (PolyPair (s, t)) = 
        PolyPair (polyadd p s, polyadd q t)

    -- (p, q) - (s, t) = (p-s, q-t)
    (PolyPair (p, q)) - (PolyPair (s, t)) = 
        PolyPair (polysub p s, polysub q t)

    -- (p, q) * (s, t) = (p*s - (1-x^2)*q*t, p*t + q*s)
    (PolyPair (p, q)) * (PolyPair (s, t)) = 
        let ps = polymult p s
            qt = polymult q t
            pt = polymult p t
            qs = polymult q s
            term2 = polymult oneMinusX2 qt
        in PolyPair (polysub ps term2, polyadd pt qs)

    -- Lifting an integer constant c to (c, 0)
    fromInteger n = PolyPair ([fromInteger n], [0])

    -- Abs and signum are typically required for Num but often left undefined 
    -- for complex-like polynomial structures.
    abs = error "abs not defined for PolyPair"
    signum = error "signum not defined for PolyPair"

-- Define the polynomial (1 - x^2) as a list of coefficients: 1 + 0x - 1x^2
oneMinusX2 :: [Double]
oneMinusX2 = [1, 0, -1]

poly1 :: PolyPair -> [Double]
poly1 (PolyPair (p, _)) = p

poly2 :: PolyPair -> [Double]
poly2 (PolyPair (_, q)) = q


--------------------------------------------------------------------------------
------------- Helpers ----------------------------------------------------------
--------------------------------------------------------------------------------

-- isEven        M.I.A
-- parityInt 

-- Typical degree calculation for dsp-style lists [a0, a1, ..., an]
polyDegree :: (Num a, Eq a) => [a] -> Int
polyDegree [] = -1  -- The zero polynomial has degree -1 by convention
polyDegree as = length (reverse (dropWhile (==0) (reverse as))) - 1


getLeadingOrderTerm :: (Eq a, Num a) => [a] -> Maybe (a,Int)
getLeadingOrderTerm coeffs = 
    case dropWhile (== 0) (reverse coeffs) of
        []    -> Nothing                -- The polynomial is zero
        (x:_) -> Just (x, length (dropWhile (== 0) (reverse coeffs)) - 1)

-- | Checks if p(x)^2 <= 1 for all x in [-1, 1] using dsp package functions
checkPolyBounds :: [Double] -> Bool
checkPolyBounds [] = True
checkPolyBounds coeffs = 
    let 
        epsilon = 1e-9
        
        -- 1. Use dsp's built-in derivative function
        pPrimeCoeffs = polyderiv coeffs

        -- 2. Find roots of the derivative to locate stationary points
        pPrimeRoots = roots epsilon 1000 (map (:+ 0) pPrimeCoeffs)

        -- 3. Filter for real roots strictly within the interval (-1, 1)
        criticalPoints = [ realPart r | r <- pPrimeRoots
                         , abs (imagPart r) < epsilon
                         , realPart r > -1 + epsilon
                         , realPart r < 1 - epsilon ]

        -- 4. Points to evaluate: endpoints and internal critical points
        testPoints = [-1.0, 1.0] ++ criticalPoints

        -- 5. Use dsp's built-in evaluation function
    in all (\x -> abs (polyeval coeffs x) <= 1.0 + epsilon) testPoints
