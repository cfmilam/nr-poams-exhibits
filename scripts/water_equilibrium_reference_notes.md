# Water: a bounded correlated equilibrium reference

16 September 2026. **Completed finite-basis local minimum; not a basis-limit or native-closure validation.**

The atomic-to-molecular reference now goes beyond a slope at a fixed geometry.
All three internal coordinates were optimized independently in a declared
effective model, and the stationary point passed full local-curvature,
half-step derivative, reference-stability and normalization checks.
The larger basis still has a nonzero gradient there. That matters: the result
is a resolved local minimum of the stated approximation, not a precision
solution for physical water.

## Result and what selected it

| Quantity | Earlier calculated starting point | Completed T-basis reference |
|---|---:|---:|
| First O–H distance | 94.13526 pm | 96.05707 pm |
| Second O–H distance | 94.13528 pm | 96.05707 pm |
| Included angle | 108.01098° | 104.28913° |

Equal distances were **not imposed**. Each distance and the included angle were
free coordinates. No observed water geometry or target energy was supplied to
the optimizer, used to tune its parameters, or used as a stopping criterion.
This is one local search from a previously calculated geometry, not a proof of
a global minimum. Printed digits identify the numerical result; they are not
experimental accuracy claims.

Relative to that starting point, the same T-basis objective decreases by
0.00092087984 E*, or **0.00605910 PHz in E/h units**. This is an energy-equivalent
rate difference between geometries—not an emitted spectral line, a binding
energy, or a physical internal clock. The complete minimum energy is
−76.39965663296607 E* with the model's fixed energy zero.

## The declared model

The calculation is all-entry restricted coupled cluster with singles, doubles
and perturbative triples, **RCCSD(T)**. It correlates all ten mobile entries in
the complete supplied molecular space: 105 spatial basis functions for T.
The H/O Gaussian basis definitions and occupied-reference exports are the
published atomic inputs. Their overlaps, normalization and occupied support
are checked again in the molecular representation.

The operator is the conventional nonrelativistic inverse-distance operator
with infinite-inertia point compact components, read through the calibrated
NR–POAMS account. Measured alpha supplies the common units. This is an adopted
numerical reference, **not an independent derivation of that operator or of
native molecular closure from the ontology**. There is no density fitting,
frozen-core approximation, empirical geometry correction, added isolated-atom
correlation energy, or appended dispersion tail. Explicit four-centre
integrals are used; this does not mean PySCF's optional integral-direct CC
algorithm was enabled.

CCSD(T) contains perturbative triples, not complete triples or quadruples.
The method definition and restricted/unrestricted distinction follow the
[primary PySCF documentation](https://pyscf.org/user/cc.html). Single-reference
diagnostics screen for evident difficulty; they do not bound missing
correlation or turn this into the explicitly correlated helium benchmark.

## The comparison that justified optimization

Before optimization, a separate registered 30-point calculation used the
**same old geometry** for the D/T/Q cardinal ladder and a T representation with
an additional oxygen diffuse extension. Function counts are 45/105/201/137.
RHF, MP2, CCSD and CCSD(T) energies were retained at each point.

| Representation | CCSD(T) angular slope, (E/h)/rad in PHz/rad | Symmetric stretch slope, E*/Å |
|---|---:|---:|
| D | +0.05534999 | not tested |
| T | +0.05258013 | −0.07027687 |
| Q | +0.05093668 | −0.06097554 |
| T + oxygen diffuse extension | +0.05303277 | not tested in this comparison |

Positive angular slopes favor a smaller angle at this particular geometry;
negative symmetric-stretch slopes favor longer bonds. All seven registered
local-descent staging checks passed. The 25% staging ceilings were permission
to conduct a **bounded local search**, not 25% accuracy estimates or a
basis-convergence certificate.

Two-step finite differences were checked against independent analytic RHF and
CCSD gradients at T/Q. Those analytic checks are **not analytic CCSD(T)
gradients**; the total CCSD(T) derivatives use energy differences. Restricted
versus unrestricted closed-shell controls, explicit pair contractions,
rigid-chart/permutation and energy-zero checks also passed. The original
comparison consumer audit records 306 checks.

## The separately registered equilibrium experiment

The optimizer started only after the comparison audit and staging gates
passed. Its source, input states, comparison and registration were frozen
before the first objective call. L-BFGS-B used coordinate scales
(0.05 Å, 0.05 Å, 0.2 rad); bounds were the old point ±0.15 Å for each distance
and ±0.35 rad for the angle. No bound was reached. The cap was 20 iterations,
200 distinct energy points and two hours, on one thread.

The search converged in five iterations and six objective/gradient evaluations.
Validation brought the distinct objective-point count to 87, plus one
independently solved final-center state. All completed by 08:09 CDT; no optimizer
remains active for this release.

- Richardson gradients are approximately (1.49e−7, 1.44e−7, 3.16e−9) in
  E*/Å, E*/Å, E*/rad, below the registered 2e−5 ceilings. Full/half-step
  disagreement stays below 7.7e−8, against the 5e−6 ceiling.
- The full internal Hessian includes cross terms. Its coordinate-scaled
  eigenvalues are (0.00440261, 0.00485353, 0.00679822) E*, all above the
  registered 1e−4 floor. Half-step diagonal changes pass independently.
  These are positive local curvatures, **not vibrational frequencies**.
- Swapping the two endpoints changes the energy by 5.4e−13 E*, without
  imposing bond equality. The final reference passes internal/external RHF
  stability screening, retained support and ten-entry normalization.
- Final T1/D1 diagnostics are 0.00875/0.02023, within the registered screening
  bounds; their passing is not an error estimate.
- The original equilibrium consumer audit records 477 checks. A separate
  publication review independently reconstructs stencils, cross curvatures,
  half-step diagonals, endpoint exchange, staging gates and bounds, rebuilds the
  final overlap and RHF/CCSD energy contractions, and recomputes the final
  perturbative triples from the saved amplitudes. All 23 frozen input hashes
  remain unchanged. This is numerical cross-checking, not a second physical
  method or a second electronic-structure implementation.

## The basis limitation is visible, not averaged away

At the T minimum the Q gradients are approximately
(+0.00451533, +0.00451533, −0.000247258) in E*/Å, E*/Å, E*/rad.
They are not zero: **the T minimum is not the Q minimum**.

Applying the T Hessian to the independently evaluated gradients gives these
local Newton displacement estimates:

| Change of representation at the T minimum | Each distance | Angle |
|---|---:|---:|
| T → Q | −0.24845 pm | +0.19203° |
| T → oxygen triple-diffuse extension | +0.00550 pm | −0.02108° |

Neither row is a relaxed geometry, an uncertainty interval, an extrapolation
to the complete basis, or a new stationary point. The Q result motivates a
separately registered Q relaxation and subsequent convergence checks; this
release does not automatically launch one. The smaller diffuse response
does not remove cardinal-basis sensitivity.

No full higher excitations, compact recoil, relativity or nuclear-motion
corrections were included. No dissociation curve, atomization energy,
spectroscopic prediction or experimental-accuracy qualification was attempted.
Oxygen's published **1.018046% atomic axial cardinal sensitivity** and the
separate atomic removal-energy limitations remain in the atomic ledger;
they are not copied onto a molecular geometry as an invented percentage error.

## Prior results are not rewritten

The native water failure, methane partial result, stopped methanol run,
failed radial conversion and failed s/p occupied-seed restriction remain.
This fuller conventional reference is not a rescue of those native claims,
nor evidence of a minimal-space computational advantage. The previous
[fixed-geometry record](molecular_reentry_notes.md) and its failure trail stay
available unchanged.

## Reproduction and custody

[Numerical summary](water_equilibrium_reference_results.json) ·
[source, registrations, inputs and scientific evidence](water_equilibrium_reference_record.zip).

The archive uses workspace-relative paths and a SHA-256 manifest. It contains
the two new registrations, frozen sources, all 30 comparison records, all 87
optimization/check points, saved center states, original audits, publication
review, needed atomic inputs, and the earlier re-entry evidence. Runtime
checkpoints, scheduler payloads, private conversations and process samples are
excluded. The frozen sources are unchanged.

Environment: Python 3.13.15, PySCF 2.14.0, NumPy 2.5.3 and SciPy 1.18.1;
one BLAS thread. The original comparison freeze accidentally labels the PySCF
package namespace as its version; the separately preserved environment receipt
records the actual version without altering the freeze.

After extracting to a **new** directory, run the original consumer audits via
the read-only replay wrapper, then the publication reviewer:

```sh
python poams-review/artifacts/water-publication-review-2026-09-16/audit_saved.py
python poams-review/artifacts/water-publication-review-2026-09-16/review.py
```

These check saved evidence; they do not run an optimizer. The wrapper executes
the original audits unchanged but redirects their audit-report writes into the
review directory, preserving the frozen original receipts. It verifies their
check counts, check names and successful outcomes. For a fresh numerical rerun, use a
second extracted copy, preserve registrations and original results elsewhere,
and select new output directories before executing the producer scripts. Do
not run the old failed s/p attempt as an automatic replacement job, and do not
overwrite the archived first-attempt failures.
