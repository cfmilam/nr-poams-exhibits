# Molecules resumed: fixed-geometry water reference, 16 September 2026

**Completed bounded diagnostic, not a new equilibrium validation.**
The atomic library now supports a directly reproducible molecular reference.
This test uses the exported H/O Gaussian basis definitions and their occupied
reference directions; it does not add isolated atomic correlation energies to
a molecular energy or claim a new derivation of the conventional operator.

## Two attempts, separately preserved

The first preregistered attempt saved five smaller-basis points, then stopped
when a proposed s/p restriction excluded small nonzero parts of oxygen's
occupied reference. The two coefficient norms outside the restricted space
were 0.0120231 and 0.00587588. These chart-dependent norms are not probabilities.
The test's failed occupied-seed gate was not weakened or relabeled a pass.

A separately registered full-space reference then retained **every** resolved
direction. It contains 15 main points, two energy-zero/chart controls, and six
atomic-fragment own/ghost-basis diagnostics. All completed. No reduced/minimal
closure or computational speed advantage is established.

## What was held fixed

The geometry is the earlier water calculation's minimum: O-H distances
0.9413525667 and 0.9413527642 angstrom, angle about 108.01097747 degrees.
It is not the experimental equilibrium geometry. Five angles per basis are
this value plus 0, ±0.5 and ±1 degree, with both distances fixed. No optimizer,
experimental scoring, fitted correction or inferred equilibrium angle enters.

The model is nonrelativistic, infinite-inertia point-compacts, with measured
alpha supplying the common length and energy units. Closed-shell RHF and
all-entry MP2 use the same direct four-centre integrals in each representation.
MP2 is a second-order joint-configuration diagnostic, **not UCCSD(T)** and not
the earlier explicitly correlated helium calculation.

## Local angular slopes

Positive means that increasing the angle raises the energy locally. Rate-first
units are derivatives of E/h in PHz per radian; divide by 6.579683920428908 for
the dimensionless E* energy slope. Values use the two-step Richardson estimate.

| Atomic representation | Functions | RHF | RHF + MP2 |
|---|---:|---:|---:|
| Single-augmented D | 45 | 0.03894678 | 0.05531365 |
| Single-augmented T | 105 | 0.03370184 | 0.05224692 |
| T, oxygen triple-diffuse extension | 137 | 0.03377708 | 0.05265192 |

Both descriptions favor a smaller angle at the tested point. Within each
fixed representation, the MP2 change in slope exceeds ten times the summed
0.5/1-degree stencil discrepancy and the registered absolute floor. This is
a resolved same-space correlation effect, not an optimized angle or proof
that correlation alone explains the earlier water error. D/T changes combine
radial and angular improvements; the failed restriction supplies no isolated
angular-only comparison. Nor do two cardinal sizes establish a basis limit.

## Checks and remaining limits

- Five imported atomic states reproduce their stored overlaps and normalization.
  Every occupied direction is retained in the full-space calculation.
- Independent unique double-replacement contractions reproduce library MP2
  energies. Pair gaps are positive; the registered perturbative diagnostics pass.
  These tests do not bound the missing higher-order correlation.
- At T, analytic minus finite-difference angular gradients are -1.8e-10 E*/rad
  (RHF) and +3.3e-8 E*/rad (MP2). Rigid rotation/translation plus endpoint
  permutation, and a known one-body energy-zero shift, change corrected energies
  by less than 7e-12 E*. Positive full overlap support is retained without an
  energy-dependent truncation; no general ill-conditioned-chart theorem follows.
- T atom-fragment own-minus-ghost basis energies sum to 0.000382831 E* at
  mean-field and 0.003055308 E* at MP2 (rate equivalents about 0.002519 and
  0.020103 PHz). These are basis-superposition diagnostics, **not bond energies**.
  Oxygen's spin sector is retained, but its ghost-induced orbital orientation
  is not independently certified; this comparison cannot qualify a dissociation
  curve or precision atomization energy.
- An independent consumer audit passed 127 custody, decimal-stencil, density,
  normalization, numerical-invariance and failure-retention checks.
- Atomic oxygen's 1.018046% axial cardinal change and the published removal-
  energy limits remain. Static atomic response is not a full excitation spectrum.
- Water's previous failed validation, methane's partial result and stopped
  methanol experiment remain unchanged. No new long geometry run is active.

The next molecular step is same-geometry higher-order correlation and further
basis checks, followed by a separately frozen equilibrium experiment only if
its numerical and approximation limits justify it.

## Reproduction

[Numerical summary](molecular_reentry_results.json) ·
[source/evidence archive](molecular_reentry_record.zip).

The archive has a workspace-relative layout and a SHA-256 manifest. It includes
both registrations, both attempts' source/logs/outputs, the five atomic JSON/NPZ
inputs, the saved starting geometry and the independent pair-contraction code.
No executable pickle arrays are required. Environment: Python 3.13, PySCF 2.14.0,
NumPy and SciPy; one BLAS thread. The run requested 3 GB and a 30-minute cap.

To reproduce without overwriting evidence, extract to a new directory, then
create a new sibling directory under `poams-review/artifacts/`, copy
`run_full.py` and `audit.py` from `molecular-reentry-full-2026-09-16/` into that
new directory, and run them there in that order. They use the immutable input
paths in the extracted tree and write only into the new sibling directory.
Preserve the original failed attempt as part of the audit, not as a restart job.

Primary implementation references: [PySCF MP2](https://pyscf.org/user/mp.html)
and [molecular structure and basis API](https://pyscf.org/user/gto.html).
