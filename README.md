# ⚛️ Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) Gauge–Higgs Unification

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21432626-1B6F8C?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.21432626)
[![License](https://img.shields.io/badge/License-Apache_2.0-B5530F)](LICENSE)
[![Lean 4](https://img.shields.io/badge/Lean_4-axiom--free-2C2C2C)](https://leanprover.github.io/)
[![Language](https://img.shields.io/badge/paper-EN_%2B_ES-1B6F8C)](paper/)

**📄 Paper (EN + ES), Lean certificates, the reusable tool & all data on Zenodo → https://doi.org/10.5281/zenodo.21432626**

> ### 📚 Part **I** of a series
> - **Part I — *Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) GHU*** (this repo): embeds the SM quarks in $(3,\mathbf{60})$ and proves the tadpole never obstructs.
> - **Part II — *Three Gates to a Quark Generation***
>   → [github.com/karlesmarin/su4-sm-cell-criterion](https://github.com/karlesmarin/su4-sm-cell-criterion) · [Zenodo 10.5281/zenodo.21432627](https://doi.org/10.5281/zenodo.21432627)
> - **Part III — *A Centre-Charge Selection Rule for the Wilson-Line Potential***
>   → [github.com/karlesmarin/centre-parity-selection](https://github.com/karlesmarin/centre-parity-selection) · [Zenodo 10.5281/zenodo.21438226](https://doi.org/10.5281/zenodo.21438226)
> - **Part IV — *Schur Functions at $(1,-1,t,t^{-1})$***
>   → [github.com/karlesmarin/schur-nonidentity-o4](https://github.com/karlesmarin/schur-nonidentity-o4) · [Zenodo 10.5281/zenodo.21463000](https://doi.org/10.5281/zenodo.21463000)

An **exact, machine-checked consistency analysis** of embedding a colored quark block in the
dimension-60 representation $(3,\mathbf{60})$ of $SU(4)$ on the orbifold $T^2/\mathbb{Z}_2$, over the
two-Higgs-doublet gauge–Higgs model of Akamatsu, Hirose, Maru and Nago (AHMN, [arXiv:2603.05857](https://arxiv.org/abs/2603.05857)).
The gauge/Higgs sector is theirs; the **fermion embedding is ours**. Physics first; every consistency
condition is a reproducible certificate, not a case-by-case argument. The paper is honest about scope:
the physical **viability verdict is left open**, and we say so throughout.

## 🧭 The journey — three questions, each stricter than the last

1. **Does $(3,\mathbf{60})$ contain the Standard Model?** — Yes. Under the unique unbroken hypercharge
   $Y=(-7q_8+4q_{15})/18$ it holds $Q_L,u_R,d_R$ with their physical charges. An exhaustive scan shows
   $(3,\mathbf{60})$ is the **minimal** $SU(4)$ irrep (with its conjugate) able to host a full quark
   generation.
2. **Is the resulting theory consistent?** — Yes, at the orbifold level. We exhibit a **45-multiplet
   completion** that cancels the eight bulk/zero-mode anomaly coefficients **and** the localized tadpole
   $\Theta_{\rm base}=(126,-36)$ **together**, with integral localized cubic residues. The localized-anomaly
   fragment is certified **axiom-free in Lean 4**.
3. **Is it physically viable?** — **Open.** Lifting the exotics without lifting $Q,u_R,d_R$ is an exact
   boundary-mass matching problem, left open pending the full brane spectrum; each outcome carries a
   concrete low-energy signature (a Kaluza–Klein/2HDM spectrum, or a light exotic / vector-like quark).

## 🔬 The main result

Contrary to the naive worry that a nonzero localized tadpole — chirality-independent, so *not* fixable by
vector-like pairs — obstructs every anomaly-free completion, **anomaly freedom and localized-tadpole
cancellation are compatible** for $(3,\mathbf{60})$, once the completion is posed and solved as an exact
integer problem. Independent verification in exact rational arithmetic gives

```
Σ (8 anomalies) = 0 ,   Σ Θ = (-126, 36) = -Θ_base ,   [SU3_C³, SU3_123³]_res = [2, 6] ∈ ℤ²
```

**And it is structural, not luck.** We prove an **Existence theorem**: on the fermion library the
localized tadpole is linearly independent of the anomalies (the combined map has rank $10=8+2$) and
freely adjustable (its image on the anomaly-neutral non-negative cone is all of $\mathbb{R}^2$), so
*every* anomaly-free completion of the block admits a tadpole-compatible one — the 45-multiplet witness
realizes a guarantee, and is representative rather than exceptional. The certificate is machine-checked
in **Lean 4** (`R60JointExistence.lean`, depending only on `propext`).

## 🧰 A reusable tool: `ghu_preflight.py`

The proof is a **rank + cone** test that is model-independent. We ship it as an exact decoupling oracle:
give it any orbifold fermion library (representations, colors, parities, chiralities) and a base block,
and it returns the exact verdict — **`UNOBSTRUCTED`** / **`NO-GO`** (with the conserved invariant φ∘Θ) /
**`CONE-OBSTRUCTED`** — together with its rank/Farkas certificate; the physics model is a single pluggable
callback. Ideal as an **exact preflight** before an expensive search (decision-complete on both sides),
and as a **landscape classifier** to hunt the genuinely obstructed models, where anomaly freedom fails to
imply tadpole freedom.

```
sage scripts/ghu_preflight.py      # runs the SU(4)/T²Z₂ (3,60) model -> VERDICT: UNOBSTRUCTED
```

## 🔒 Machine-checked in Lean 4 (axiom-free)

Two certificates, both **`sorry`-free** and depending **only on `propext`** — `#print axioms` shows no
`Classical.choice`, no `sorryAx`. Each is an exact integer identity discharged by `decide`, so the check
is the theorem, not a numerical approximation.

- **`lean/R60LocalizedAnomaly.lean`** — the *localized* (fixed-point) content of the $(3,\mathbf{60})$
  block: the localized non-Abelian cubic residues $[SU3_C^3, SU3_{123}^3]=[2,6]$, the dimension count,
  and the Wilson-line (in)dependence of the pure and mixed channels.
- **`lean/R60JointExistence.lean`** — the **Existence theorem**'s engine: on the fermion library the
  combined (anomaly, tadpole) map has rank $10=8+2$ (the tadpole is *independent* of the anomalies) and
  its image on the anomaly-neutral non-negative cone is all of $\mathbb{R}^2$ (four Farkas witnesses),
  so **every** anomaly-free completion is tadpole-compatible — the 45-multiplet witness realizes a
  guarantee, not a coincidence.

Both build in the fast `godsil-gutman-lean` environment; the axiom audit is part of each certificate.

## ✅ Reproducibility

Every number regenerates from the ancillary artifacts:

- `lean/R60LocalizedAnomaly.lean` — axiom-free, `sorry`-free Lean 4 certificate of the localized cubic
  residues, the dimension count, and the Wilson-line (in)dependence of the pure/mixed channels.
- `scripts/verify_formulas.sage` — the full formula audit (SM charges, the 8 coefficients, the witness
  cancellation) in exact rational arithmetic.
- `scripts/minimality_scan.sage` — proof that $(3,\mathbf{60})$ and its conjugate are the only $SU(4)$
  irreps of dimension ≤ 63 hosting a full quark generation.
- `scripts/r60_verify_joint.sage` — independent re-verification of the 45-multiplet witness (physical
  hypercharge convention).
- `scripts/r60_cpsat.sage` — the OR-Tools CP-SAT joint completion search.
- `supplement/` — the full 45-multiplet spectrum (`witness_45_spectrum.csv`), the raw audit log
  (`formula_audit.log`), and a reproducibility index (`SUPPLEMENT.md`).

## 📁 Contents

```
paper/    ghu_completion.pdf/.tex (EN) · ghu_completion_es.pdf/.tex (ES) · figs/
lean/     R60LocalizedAnomaly.lean · R60JointExistence.lean   (both axiom-free, propext only)
scripts/  ghu_preflight.py · structural_theorem.sage · verify_formulas · minimality_scan · r60_verify_joint · r60_cpsat · make_paper_graphics(+_es)
supplement/  SUPPLEMENT.md · witness_45_spectrum.csv · formula_audit.log
```

## 📜 License & citation

Released under [Apache 2.0](LICENSE). If you use this, please cite the Zenodo record:
**[10.5281/zenodo.21432626](https://doi.org/10.5281/zenodo.21432626)**. Author: **Carles Marín**
(independent researcher). The exact computations were carried out and cross-checked with Claude
(Anthropic) as an AI research assistant against a common machine-verifiable ground truth (Lean, exact
algebra).
