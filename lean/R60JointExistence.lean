/-
Author: Carles Marín (Claude, Anthropic, as AI research assistant).

R60JointExistence.lean — an axiom-free certificate for the Existence theorem of the
six-dimensional SU(4) gauge–Higgs (3,60) paper: the localized tadpole never obstructs an
anomaly-free completion of the block.

Each "species" of added matter contributes an exact integer 10-vector
  [ 108·(eight anomaly coefficients) , Θ01 , Θ23 ]
(the anomalies are scaled by their common denominator 108 to clear fractions; SageMath).
We exhibit four anomaly-neutral non-negative additions — sums of species — that realize
the four signed tadpole axes ±Θ01, ±Θ23 (Farkas certificates that 0 is interior to the
achievable tadpole cone), and check that two of their tadpoles are linearly independent
(nonzero 2×2 determinant), so Θ has rank 2 on ker A and no invariant φ∘Θ is conserved.
Everything is decided by the kernel: no Mathlib, no `sorry`, no extra axioms.

Build (fast, no Mathlib):  lake env lean R60JointExistence.lean
-/

namespace R60JointExistence

abbrev Vec := List Int   -- length 10: [a0,…,a7, Θ01, Θ23]

def add (u v : Vec) : Vec := List.zipWith (· + ·) u v
def sumv (l : List Vec) : Vec := l.foldl add (List.replicate 10 (0 : Int))

-- the twelve species making up the four cone witnesses (108·anomalies | tadpole)
def s86  : Vec := [-3780, -88128, -11340,  648, -288,  108, -21828, -1728, 144,  0]
def s87  : Vec := [ 3780,  88128,  11340, -648,  288, -108,  21828,  1728, 144,  0]
def s82  : Vec := [-3780, -88128, -11340, -324,  288, 1188,  11028,  1728,-144,  0]
def s83  : Vec := [ 3780,  88128,  11340,  324, -288,-1188, -11028, -1728,-144,  0]
def s129 : Vec := [ 2160,  -3564,   6480,  108, -126, -621,  -2253,  -756, -54, 12]
def s133 : Vec := [ 2160,  -3564,   6480, -108,   90,  378,    366,   540,  54, 12]
def s234 : Vec := [-2160,   3564,  -6480,  108,  -90, -378,   -366,  -540,  54, 12]
def s238 : Vec := [-2160,   3564,  -6480, -108,  126,  621,   2253,   756, -54, 12]
def s122 : Vec := [    0,   1188,  -2160,    0,    0,  126,    122,   180, -18, -4]
def s127 : Vec := [    0,  -1188,   2160,    0,    0,  207,    751,   252,  18, -4]
def s216 : Vec := [    0,   1188,  -2160,    0,    0, -207,   -751,  -252,  18, -4]
def s221 : Vec := [    0,  -1188,   2160,    0,    0, -126,   -122,  -180, -18, -4]

-- the four Farkas witnesses (each a non-negative sum of species)
def wP01 : Vec := sumv [s86, s87]              -- +Θ01
def wM01 : Vec := sumv [s82, s83]              -- −Θ01
def wP23 : Vec := sumv [s129, s133, s234, s238] -- +Θ23
def wM23 : Vec := sumv [s122, s127, s216, s221] -- −Θ23

def anomZero (v : Vec) : Bool := (v.take 8).all (· == 0)   -- all eight anomalies cancel
def tad (v : Vec) : Int × Int := (v.getD 8 0, v.getD 9 0)  -- the tadpole (Θ01, Θ23)

/-- **Existence certificate.** The four additions are anomaly-neutral and realize the four
signed tadpole axes; two of their tadpoles, (288,0) and (0,48), have nonzero determinant,
so Θ has rank 2 on the anomaly kernel and no invariant φ∘Θ is conserved. Hence every
anomaly-free completion of (3,60) admits a tadpole-compatible one. -/
theorem existence_certificate :
    (anomZero wP01 = true ∧ anomZero wM01 = true ∧
     anomZero wP23 = true ∧ anomZero wM23 = true)
  ∧ (tad wP01 = (288, 0)  ∧ tad wM01 = (-288, 0) ∧
     tad wP23 = (0, 48)   ∧ tad wM23 = (0, -16))
  ∧ ((288 : Int) * 48 - 0 * 0 ≠ 0) := by decide

#print axioms existence_certificate

end R60JointExistence
