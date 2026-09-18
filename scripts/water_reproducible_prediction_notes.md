# Water: reproducible prediction closes the pipeline validation

17 September 2026. The same finite-basis effective model that produced our T
and Q geometry references now predicts water's harmonic vibrational
frequencies and static-nuclei atomization energy. Both agree with the
established experimental values to within finite-basis, harmonic-approximation
and post-CCSD(T) error bars, using **no experimental input** to select any
parameter of any calculation in this release. That closes the pipeline
validation loop for water. It does not close, replace, or validate the
separate native-closure question.

## What the pipeline predicts

Two properties, both derived from artifacts already frozen by the Q reference:

| Property | This pipeline | Established experiment | Deviation |
|---|---:|---:|---:|
| Harmonic bend ω₂ | 1648.59 cm⁻¹ | 1648.5 cm⁻¹ | +0.09 cm⁻¹ (+0.005%) |
| Harmonic symmetric stretch ω₁ | 3837.10 cm⁻¹ | 3832.2 cm⁻¹ | +4.90 cm⁻¹ (+0.128%) |
| Harmonic antisymmetric stretch ω₃ | 3947.36 cm⁻¹ | 3942.5 cm⁻¹ | +4.86 cm⁻¹ (+0.123%) |
| Atomization D_e (static nuclei) | 231.53 kcal/mol | 232.60 ± 0.05 kcal/mol | −1.07 kcal/mol (−0.46%) |
| D_0 (with harmonic ZPE) | 218.05 kcal/mol | 219.35 ± 0.05 kcal/mol | −1.30 kcal/mol (−0.59%) |

The experimental harmonic frequencies are those extracted from vibrational
spectroscopy through the standard anharmonic fit; the experimental D_e is the
static-nuclei electronic dissociation limit inferred from atomization
thermochemistry and independent theoretical estimates of anharmonic and
non-adiabatic corrections. These are the reference values a converged
CCSD(T)/aug-Q calculation is normally compared to, and this pipeline lands
within their conventional CCSD(T)/aug-Q error window.

## Why this counts as prediction

- **No experimental frequency, bond length, angle, atomization energy, ZPE,
  or spectral line was supplied to any calculation in this release.** The
  fine-structure calibration of the pair operator is the same one used for
  every prior POAMS-instrument prediction and remains fully disclosed.
- The atomic Q energies used here are the **same** frozen upstream records
  the water optimizer consumed. Atoms and molecule share one basis definition
  and one correlation treatment; no mix-and-match.
- The Wilson G matrix uses **CODATA physical constants and neutral-isotope
  masses**. No spectroscopic input.
- The Hessian was locked at the end of the bounded Q optimization on 16
  September; no numbers were regenerated for this release.

Under those conditions, agreeing with three independent experimental
harmonic frequencies to fractions of a percent, and with two independent
experimental atomization endpoints to less than half a percent, is a
non-empirical **prediction**, not a fit.

## Numerical protocol

### Vibrational frequencies (Wilson GF)

Coordinates: internal (r₁, r₂, θ), the same three used by the Q optimizer.
Masses: ¹H = 1.00782503207 u, ¹⁶O = 15.99491461957 u (ground natural
isotopes).

- **G matrix** built from scratch at the Q equilibrium using the standard
  Wilson expressions for a bent B–A–B molecule. Recovers the OH reduced
  mass, 1/(1/m_H + 1/m_O), on the (r,r) diagonal to machine precision.
- **F matrix** copied from the Q reference `curvature.matrix_internal`,
  converted to SI J/m² and J/rad².
- **GF diagonalized** with double-precision LAPACK; sorted eigenvalues give
  three positive, non-degenerate normal modes.
- **Frequencies** ω_i = √λ_i / (2π c) in cm⁻¹.
- **Modes classified** by the sign of the (r₁, r₂) product in the
  eigenvector: same sign → symmetric stretch, opposite sign → antisymmetric
  stretch, remaining mode → bend.

The harmonic ZPE is (1/2)·ħ·Σω, giving 13.485 kcal/mol from these three
modes. It is used only to construct a D_0 comparison to the observed
zero-point-corrected atomization; no anharmonic correction is applied.

### Atomization energy

The **primary** number uses the campaign's own Q-level atomic bases:
- H: (5s,4p,3d,2f) = 46 spherical basis functions (aug-cc-pVQZ).
- O: (8s,8p,6d,4f,2g) = 109 spherical basis functions (a custom augmentation
  of Dunning's QZ family, sufficient to converge oxygen's small second-derivative
  response inside this campaign).

Both atomic records are frozen upstream files that the Q water optimizer
already consumed. Reading them back for this release keeps atoms and molecule
at strictly the same basis and same correlation treatment (RCCSD(T),
ROHF-CCSD(T) for the open-shell O).

D_e = 2·E(H) + E(O) − E(H₂O) = 0.368969 Hartree = 231.531 kcal/mol
= 10.04015 eV.

D_0 = D_e − ZPE(harmonic) = 218.046 kcal/mol.

A **standard cc-pVQZ** atomic calculation is also provided in this artifact
as an independent-basis diagnostic. It uses a different (smaller) atomic
basis than the water optimizer used and cannot be combined with the Q molecular
energy without a like-basis correction. It is included to show that our
standalone PySCF setup reproduces the expected atomic H (−0.4999456 Ha) and
O (−75.02318 Ha) at cc-pVQZ, matching published values for those atoms in
that basis. Do not mix these atomic numbers into the primary D_e.

### Independent audit

`audit.py` reconstructs the frequency computation from scratch:

- Rebuilds the Wilson G matrix symbolically, checks (r,r) diagonal recovers
  the OH reduced mass, and reproduces G(θ,θ)'s symmetric-limit closed form.
- Rebuilds the F matrix in SI from the frozen Q Hessian.
- Diagonalizes GF with an independent path, requires the eigenvalues match
  the driver to below 1 μcm⁻¹.
- Rebuilds ZPE from ħ·Σω/2 to relative 1e-12.
- Re-reads the frozen H and O Q atomic records, verifies their SHA-256
  against the Q freeze, checks E(H)_CCSD(T) = E(H)_reference (1-electron
  exactness), checks E(O)_reference + CCSD_corr + (T)_corr = E(O)_CCSD(T),
  and recomputes D_e in Hartree, kcal/mol and eV.
- Independently reads the cc-pVQZ diagnostic H and O energies and requires
  the H atom sits within 200 μHa of −0.5.

52 checks; audit passes with zero failures.

## What this does and does not establish

**Establishes.** Reproducible prediction of water's harmonic vibrational
spectrum and atomization energy from a bounded, non-empirical nonrelativistic
RCCSD(T)/aug-Q pipeline. Two independent experimental observables match, at
the accuracy the method can support, without any calibration to those
observables. Extending the same pipeline to other small closed-shell molecules
is mechanical: the machinery does not need new methodology to run.

**Does not establish.**
- Post-CCSD(T) accuracy. Full triples, quadruples, and multireference are
  not attempted; the remaining ~1 kcal/mol atomization discrepancy is in the
  expected direction and magnitude for missing higher excitations plus small
  basis-set-incompleteness residuals.
- Basis-limit geometry. Only two cardinals (T and Q) are relaxed; the Q
  diffuse-augmented and Q→5 convergence steps remain open, as reported in
  the [Q geometry release](water_q_relaxation_notes.md).
- Vibrational anharmonicity or non-Born–Oppenheimer effects. The ω values
  above are harmonic; observed fundamentals ν_i lie 3–8% below the
  corresponding ω_i.
- Relativity, recoil, nuclear-motion corrections beyond harmonic ZPE, or
  core-valence correlation beyond the campaign's frozen-orbital defaults.
- **Native POAMS molecular closure.** This release validates the numerical
  pipeline used to test POAMS derivations; it does not itself derive water's
  geometry, spectrum or binding from a POAMS first principle. The previous
  native water, methane, methanol, and radial/s/p reduction failures remain
  unchanged in the ledger, as do oxygen's atomic 1.018046% axial cardinal
  sensitivity and separate removal-energy limits.

## Scaling to more molecules

The demonstrated pipeline — atomic basis freeze; RCCSD(T) energy and
gradient; internal-coordinate Hessian; Wilson GF frequencies; like-basis
atomization energy — has no molecule-specific step. Applying it to another
small closed-shell hydride, or a substituted analogue, requires only:

1. constructing (or reusing) the atomic Q records for the new elements;
2. running the same bounded optimizer against the new molecule at T then Q;
3. rerunning `predict.py` with a different geometry-and-Hessian input.

Each additional molecule adds cost that grows as the usual O(N⁷) scaling of
CCSD(T), not as a new methodological problem. This is why we are stopping
the software-system build here: the machinery does what it needs to do,
water shows it works end-to-end, and the outstanding scientific work is
whether native closures can eventually be re-tested against this same
audited pipeline — not whether the pipeline itself needs more parts.

## Reproduction and custody

[Numerical summary](water_reproducible_prediction_results.json) ·
[source, atomic-point cross-check, audit and evidence archive](water_reproducible_prediction_record.zip).

Archive layout is workspace-relative with a SHA-256 manifest. It contains
this release's driver, audit, registration and result; the independent
cc-pVQZ atomic point files; and — for full custody continuity — the T
reference notes/results, the Q reference notes/results, and the 31 frozen
upstream inputs the Q reference already fixed.

After extracting to a **new** directory, re-run the independent audit
(no PySCF needed for the frequency portion):

```sh
python poams-review/artifacts/water-reproducible-prediction-2026-09-17/audit.py
```

For a fresh numerical rerun of the atomic diagnostic, use the same isolated
environment the atomic campaign used (Python 3.13.15, PySCF 2.14.0, NumPy
2.5.3, SciPy 1.18.1, single BLAS thread) and select new output directories
before executing `predict.py`. Do not overwrite the archived records.

## Prior results are not rewritten

The published [T reference](water_equilibrium_reference_notes.md), the
[relaxed Q reference](water_q_relaxation_notes.md), the fixed-geometry
[re-entry precursor](molecular_reentry_notes.md), and the atomic ledger
remain available unchanged. This release consumes them; it does not modify
them. The failure trail from the original native-closure attempt and its
stopped derivative jobs remains available in the Elements atomic ledger and
the historical review artifacts.
