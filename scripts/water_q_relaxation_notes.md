# Water: relaxed Q-basis reference and T→Q cross-check

17 September 2026. **Completed larger-basis local minimum; not a basis-limit, native-closure or experimental validation.**

The T-basis reference published on 16 September was authorized against a
separately registered Q relaxation and convergence assessment. That
optimization has now finished. Two independent finite-basis minima now exist
at the same declared model and coordinates. Their agreement across cardinal
representations tightens the finite-basis picture; it does not by itself
establish a basis-limit geometry, a diffuse-function limit, or experimental
accuracy.

## Result and its comparison to T

| Quantity | Published T reference | Relaxed Q reference |
|---|---:|---:|
| First O–H distance | 96.05707 pm | **95.80968 pm** |
| Second O–H distance | 96.05707 pm | **95.80968 pm** |
| Included angle | 104.28913° | **104.48054°** |

Equal distances were not imposed. Each distance and the included angle were
free coordinates, started from the exact saved T minimum without any observed
target. The relaxed Q geometry is different from T's geometry, and different
in the direction the T Hessian predicted:

| Change of representation at the T minimum | Each distance | Angle |
|---|---:|---:|
| **T Newton-step estimate** (from T Hessian and Q gradient) | −0.24845 pm | +0.19203° |
| **Q relaxation, actual shift** | **−0.24739 pm** | **+0.19140°** |
| Actual minus estimate | 0.00106 pm | −0.00063° |

Two calculations at different basis sizes, using different tools (T's Hessian
solved analytically at T; Q's optimum found by iterative descent in a larger
representation), agree on the shift to about **one part in two hundred** on
each coordinate. That is a genuine cross-check on the T Hessian's local
predictive power inside this model, not an accidental coincidence.

The Q reference energy is **−76.42463692196 E*** with the model's fixed energy
zero. The Q relaxation from T's geometry to Q's own minimum lowers the Q
energy by 1.1617×10⁻⁵ E* (−7.643×10⁻⁵ PHz in E/h units); this is a small
model-internal relaxation, not a spectral prediction.

## The declared model

Same all-entry restricted coupled cluster with singles, doubles and
perturbative triples, **RCCSD(T)**, as the published T reference. Now with the
full **Q molecular basis: 201 spatial functions** (up from 105 at T).
Nonrelativistic inverse-distance operator, infinite-inertia point compact
components, measured alpha calibration, no density fitting, no frozen-core
approximation, no empirical geometry input, no added isolated-atom
correlation energy or dispersion tail, explicit four-centre integrals. Method
definition and restricted convention follow the
[primary PySCF documentation](https://pyscf.org/user/cc.html).

CCSD(T) still means perturbative triples, not complete triples or quadruples.
Single-reference diagnostics still screen for evident difficulty; they still
do not bound missing correlation.

## The separately registered Q experiment

The optimizer started only after the published T reference was audited and
its consumer review passed. Its source, all frozen inputs, the T result and
audit, and a new registration were hashed and frozen before the first
objective call, then re-verified before writing the final result. L-BFGS-B
used coordinate scales (0.05 Å, 0.05 Å, 0.2 rad); bounds were the T minimum
±0.15 Å for each distance and ±0.35 rad for the angle. No bound was reached.
The cap was 20 iterations, 200 distinct energy points and two hours on one
thread. The search converged in **4 iterations and 5 objective/gradient
evaluations**. Validation brought the distinct objective-point count to **80**,
plus one independently solved final-center state. All completed by
11:29 CDT; no optimizer remains active for this release.

- Richardson gradients at the Q minimum are approximately
  (2.80×10⁻⁷, 3.26×10⁻⁷, 1.66×10⁻⁸) in E*/Å, E*/Å, E*/rad, well below the
  registered 2×10⁻⁵ ceiling. Full/half-step disagreement stays below
  7.04×10⁻⁸ (ceiling 5×10⁻⁶).
- The full internal Hessian includes cross terms. Its coordinate-scaled
  eigenvalues are (0.00445, 0.00492, 0.00679) E*, all above the registered
  10⁻⁴ floor. Half-step diagonal changes pass independently. These are
  positive local curvatures, **not vibrational frequencies**.
- Independent analytic RHF and CCSD gradients at the Q minimum agree with the
  finite-difference stencils to about 4×10⁻⁹ and 1×10⁻⁷ respectively. These
  are analytic RHF/CCSD gradients used as cross-checks, **not analytic
  CCSD(T) gradients**; the total CCSD(T) derivatives use energy differences.
- Swapping the two endpoints changes the energy by −2.6×10⁻¹³ E* without
  imposing bond equality. The final reference passes internal/external RHF
  stability screening, retained support and ten-entry normalization.
- Final T1/D1 diagnostics are within the registered screening bounds; their
  passing is not an error estimate.
- The Q consumer audit records **632 checks**. It reconstructs stencils, cross
  curvatures, half-step diagonals, endpoint exchange, staging gates and
  bounds; rebuilds the final overlap and RHF/CCSD energy contractions;
  recomputes the final perturbative triples from the saved amplitudes; and
  verifies all 31 frozen input hashes remain unchanged. This is numerical
  cross-checking, not a second physical method or a second electronic-structure
  implementation.

## Diffuse response is present, small, and unresolved

At the relaxed Q geometry, adding oxygen's third diffuse extension to a T
representation changes the CCSD(T) gradient by
(−8.29×10⁻⁵, −8.29×10⁻⁵, +5.37×10⁻⁵) in E*/Å, E*/Å, E*/rad. Solving with the
Q Hessian gives a local Newton displacement of about **+5.5×10⁻³ pm** on each
bond and **−0.021°** on the angle. This is roughly **fifty times smaller** than
the T→Q shift and points in the same physical direction.

That diffuse response is a lower-cardinal comparison at a fixed geometry, not
a Q-basis diffuse-convergence proof or a complete-basis geometry. A future
release would need at least one relaxation in a diffuse-augmented Q basis
and, ideally, a Q→5 comparison before the geometry could be quoted as
basis-converged. This release does not automatically launch either.

## What still isn't established

- **Complete-basis geometry.** Q lies below T, but two relaxed cardinal
  points do not fit a rigorous extrapolation curve. The T→Q shift and the
  smaller diffuse response together **make CBS convergence plausible in
  direction**, not proven in magnitude.
- **Higher excitations.** Perturbative triples remain the ceiling. Full
  triples, quadruples, and multi-reference corrections were not attempted
  and would move the geometry in an unquantified direction.
- **Physical corrections.** No compact-recoil, relativity, nuclear-motion or
  dissociation prediction was attempted. This is a Born–Oppenheimer, static,
  point-charge electronic minimum in the adopted operator.
- **Comparison to experiment.** No observed geometry was supplied to select
  parameters, and none is claimed to be validated. Printed digits identify
  the numerical result; they are not experimental accuracy statements.
- **Native closure.** The RCCSD(T) reference is an adopted numerical
  standard against which the earlier native-closure attempt is expected to be
  eventually re-tested. It is **not itself** a rescue of that attempt, and does
  not remove the previous native water, methane, methanol or radial/s/p
  reduction failures. Oxygen's atomic **1.018046% axial cardinal
  sensitivity** and the separate atomic removal-energy limits remain in the
  atomic ledger.

## Prior results are not rewritten

The published [T reference](water_equilibrium_reference_notes.md) and its
[fixed-geometry re-entry precursor](molecular_reentry_notes.md) remain
available unchanged. The failure trail from the original native-closure
attempt and its stopped derivative jobs remains available in the Elements
atomic ledger and the historical review artifacts. This release adds a
larger-basis reference; it does not replace, retract or re-open any of them.

## Reproduction and custody

[Numerical summary](water_q_relaxation_results.json) ·
[source, registrations, inputs and scientific evidence](water_q_relaxation_record.zip).

The archive uses workspace-relative paths and a SHA-256 manifest. It contains
this Q run's source, freeze, per-point evidence, final-center state, original
audit and publication review, plus the 31 frozen upstream inputs enumerated
in the Q freeze (T reference source and result, T comparison campaign and its
outputs, the atomic H/O correlated records at cardinals 2/3/4 and third-diffuse,
and the two registered specifications). Runtime checkpoints, scheduler
payloads, private conversations, and process samples are excluded. The frozen
sources are unchanged.

Environment: Python 3.13.15, PySCF 2.14.0, NumPy 2.5.3 and SciPy 1.18.1; one
BLAS thread.

The Q reference builds on the T reference and consumes the same audit machinery.
After extracting to a **new** directory, verify custody and re-run the
independent audit:

```sh
python poams-review/artifacts/water-q-relaxation-2026-09-16/audit.py
```

That checks saved evidence; it does not run an optimizer. For a fresh
numerical rerun, use a second extracted copy, preserve registrations and
original results elsewhere, and select new output directories before executing
the producer scripts. Do not overwrite the archived first-attempt records.
