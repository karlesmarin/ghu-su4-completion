# Supplementary material

**Anomaly- and Tadpole-Compatible Fermion Completion of Six-Dimensional SU(4)
Gauge–Higgs Unification** — Carles Marín (Claude, Anthropic, as AI research assistant).

Everything reported in the paper is a finite exact computation. This supplement collects
the complete data and the raw verification logs so every number can be regenerated.

## Files

| File | Contents |
|------|----------|
| `witness_45_spectrum.csv` | The complete 45-multiplet completion (17 species) as machine-readable data: index, SU(4) Dynkin `(a,b,c)`, color rep, parities `(p0,p2)`, 6D chirality `chi6`, multiplicity. |
| `formula_audit.log` | Raw output of `verify_formulas.sage`: SM charges under the physical hypercharge, the 8 base anomaly coefficients in both the `/9` and physical `/18` normalizations, the 45-multiplet cancellation re-verified in both, tadpole and local-cubic residues, and the Q-matrix rank. |

## Reproduction scripts (repository root)

| Script | Regenerates |
|--------|-------------|
| `verify_formulas.sage` | The full formula audit (`formula_audit.log`). Exact rational arithmetic, SageMath. |
| `minimality_scan.sage` | Proof that `(3,60)=(0,2,1)` and its conjugate `(1,2,0)` are the ONLY SU(4) irreps of dimension ≤ 63 whose T²/Z₂ zero modes host a full SM quark generation under a single hypercharge — i.e. `(3,60)` is minimal. |
| `r60_verify_joint.sage` | Independent re-verification of the 45-multiplet witness (8 anomalies = 0, tadpole = 0, local cubic residue ∈ ℤ²). |
| `r60_cpsat.sage` | The OR-Tools CP-SAT joint completion search. |
| `structural_theorem.sage` | Proves `rank[A;Θ] = 10` and that the anomaly-neutral tadpole cone is all of ℝ² — the **Existence theorem**: the localized tadpole never obstructs an anomaly-free completion. |
| `ghu_preflight.py` | The reusable **rank + cone** obstruction oracle: give it any orbifold fermion library + base block, it returns `UNOBSTRUCTED` / `NO-GO` / `CONE-OBSTRUCTED` with its certificate. |
| `R60LocalizedAnomaly.lean` | Axiom-free, `sorry`-free Lean 4 certificate of the localized cubic residues, the dimension count, and the Wilson-line (in)dependence of the mixed/pure channels. |
| `R60JointExistence.lean` | Lean 4 certificate of the Existence theorem (four Farkas cone witnesses + independence), kernel-checked, depending only on `propext`. |

## The 45-multiplet completion (full list)

| # | Dynkin (a,b,c) | color | (p0,p2) | chi6 | mult |
|---|----------------|-------|---------|------|------|
| 1 | (0,0,1) | 1 | (−,−) | + | 8 |
| 2 | (0,0,3) | 1 | (+,−) | − | 3 |
| 3 | (0,0,3) | 1 | (−,+) | − | 3 |
| 4 | (0,0,4) | 1 | (−,+) | + | 2 |
| 5 | (0,1,0) | 1 | (+,+) | − | 2 |
| 6 | (0,1,1) | 1 | (+,+) | − | 2 |
| 7 | (0,1,1) | 1 | (+,+) | + | 2 |
| 8 | (1,0,1) | 1 | (+,+) | + | 2 |
| 9 | (0,0,1) | 3 | (−,+) | + | 1 |
| 10 | (0,0,1) | 3 | (−,−) | − | 4 |
| 11 | (0,0,3) | 3 | (−,−) | + | 2 |
| 12 | (0,2,0) | 3 | (−,+) | − | 2 |
| 13 | (0,0,1) | 3̄ | (+,+) | + | 5 |
| 14 | (0,0,1) | 3̄ | (−,+) | + | 3 |
| 15 | (0,0,3) | 3̄ | (+,−) | − | 1 |
| 16 | (0,0,3) | 3̄ | (−,+) | − | 1 |
| 17 | (0,1,1) | 3̄ | (−,−) | − | 2 |

**Total: 45 multiplets** (color singlets 24, triplets 9, antitriplets 12).

## Verified results (from `formula_audit.log`)

- SM quark charges under `Y=(-7q8+4q15)/18`: `Q_L=1/6`, `u_R=2/3`, `d_R=-1/3`.
- Base 8 anomaly coefficients (physical Y): `[60, -579, 180, -3, 17/6, 23/2, 2293/18, 17]`.
- 45-multiplet cancellation: 8 anomalies `= 0` (both normalizations), tadpole `Σ+base = (0,0)`,
  local cubic residue `= [2,6] ∈ ℤ²`.
- Base tadpole `Θ_base = (126,-36)`; local cubics `[-3/2, -3/4]`; Q-matrix rank 2 (det −19449).
