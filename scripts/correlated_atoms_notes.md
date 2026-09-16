# Correlated C/N/O atomic preparation

This layer upgrades the H/C/N/O molecular **controls**, not the entire
118-element library. The existing mean-field references and the exact H /
correlated finite-inertia He-4 benchmarks remain unchanged. The new Gaussian
H/He calculations are implementation and basis tests, not replacements for
those better benchmarks. Molecules remains on hold.

## What is calculated

In the Atomic Ledger's calibrated units,

    a* = hbar/(m c alpha), E* = m c^2 alpha^2,
    H/E* = sum_i [-laplacian_i/2 - Z/r_i] + sum_(i<j) 1/r_ij.

The compact component is a fixed, infinitely massive point. The supplied
antisymmetric state rule, kinetic operator and pair potential are conventional
nonrelativistic atomic mathematics, not newly proved NR axioms. Alpha and
mass calibrate units; atomic energy, size and response are not fitted.
No isotope's internal compact geometry or regional energy allocation is derived.

State inputs are C 3P, N 4S, O 3P and the removal states C+ 2P, N+ 3P, O+ 4S.
Use a spin-pure restricted-open-shell reference (ROHF), followed by separate
occupied/virtual Fock diagonalizations for each spin and symmetry block.
PySCF's unrestricted CCSD equations and perturbative triples, CCSD(T), are
then solved. The reference is ROHF; the correlation implementation is UCCSD,
not a claimed spin-adapted coupled-cluster implementation. All occupied
entries, including the inner pair, are correlated. Four-centre integrals are
evaluated directly, without the previous RI-JK factorization.

The full independent atom is solved in each basis, charge state, angular
branch and applied perturbation. Each result is checkpointed. This does not
establish a global search over every multiplet or an exact correlated solution.

## Numerical representation and controls

For C/N/O use cc-pCVDZ, cc-pCVTZ and cc-pCVQZ with the matching diffuse shells
from aug-cc-pVnZ. H/He use aug-cc-pVnZ. The extended diffuse check appends one
more shell per angular degree with exponent one third of the corresponding
diffuse exponent. Exact arrays are exported. Basis functions are imported
numerical representation data, not geometrically derived physical constants.

The completed two-layer oxygen calculation changed its response by 1.43%
relative to one layer. The completed third-layer check adds exponents at one
ninth of the original diffuse exponent, at both T and Q cardinal sizes.
Its average response changes by only 0.00594% from the second layer. This
resolves the observed diffuse-tail sensitivity, not every basis limitation:
oxygen's parallel branch still changes by 1.018% between T and Q.
Both successive cardinal and diffuse response changes must be below the
initial 1% gate for each branch as well as the average. An average cannot
clear a failing branch; matching a published response is not a substitute
for convergence. The previous solver source is archived by its record-linked
SHA-256, and the originally used one/two-layer basis arrays are unchanged.

SCF convergence is requested at 1e-12 E*, gradient 1e-8, with a final gate
below 1e-7. CC convergence is requested at 1e-11 E* and 1e-9 amplitude change;
an additional update explicitly checks the final amplitude residual below
1e-8. Spin counts, reference S^2, symmetry occupations, overlap conditioning,
orbital metric, density normalization and radial integration are checked.

H has exact nonrelativistic controls E=-0.5 E*, mean r=1.5 a*, rms r=sqrt(3) a*,
and static response 4.5. Finite Gaussian representations are compared with
these controls; their discretization error is not reclassified as new physics.
For He, CCSD agrees with independent full configuration interaction in the
same finite basis. Both are compared with the existing infinite-inertia
Hylleraas result, not the physically different finite-inertia He-4 result.

For each C/N/O neutral and removal state a separate all-entry **full
configuration interaction** calculation in cc-pVDZ checks correlation
truncation. This smaller-space difference is not transplanted to the larger
basis as a correction or represented as a rigorous larger-basis error bound.
The supplied scalar T1 norms are diagnostics, not certificates that all
multireference effects are absent.

## Response: no angular branch swapping

Use H(F)=H(0)+F sum_i z_i with fixed C2v occupations. For C and O, separately
follow the two branches corresponding to the supplied P-state axis parallel
and perpendicular to the perturbation. At zero perturbation their calculated
energies must agree. For the S-state nitrogen only one branch is needed.

    chi(F) = -[E(+F)+E(-F)-2E(0)]/F^2,
    chi(0) estimate = [4 chi(F/2)-chi(F)]/3,
    chi_average = (chi_parallel+2 chi_perpendicular)/3.

F=0.001 and 0.0005 in E*/a* units. Both signs are actually solved. The
energy's evenness and the two curvature estimates are checked; the field-step
difference must be below 0.2%. Extrapolation here removes the leading
finite-perturbation term, not basis or correlation error. These are static
angular-term responses, not J-resolved spin-orbit tensors, dynamic spectra,
or a new C6 prediction. Multiplying chi by a*^2/E* gives a mechanical
compliance; in conventional electric-field units an extra squared charge
conversion is required. Do not confuse chi with the fine-structure alpha.

The two independent branch energies are not two physical species. P states
are degenerate before the perturbation. The orientation average and axis
convention must travel with any use of the response tensor.

## Size and downloadable state convention

The size calculation uses the **CCSD left/right (lambda) one-body density**.
This is an unrelaxed-orbital density expectation, not a full CCSD(T) response
derivative or a normalized squared one-entry wavefunction. It integrates to
the total mobile-entry count. The orientation-averaged radial profile is
divided by that count, hence integrates to one. Mean, rms and 90% distances
are distinct statistics and none is a rigid surface, core radius or covalent
bond radius. The analytic r^2 integral and nested radial grid check the profile.

The main archive contains JSON for every field/state and the reproduction
scripts. Separate element archives contain NPZ for zero-field states and
their matching JSON records. Oxygen's larger package is split by diffuse
layer to keep every download below the hosting limit. The display JSON lists
all packages in `state_archives`; `state_archive` identifies the package
containing the selected result. The panel links every archived basis package.
Load NPZ using `allow_pickle=False`. `mo_coeff[spin,mu,p]` expands normalized
orbitals in the exported PySCF real-spherical, normalized contracted AO basis;
`C.T @ overlap @ C = I`. `mo_occ` is 0/1 per spin. `t1_0/t1_1` are alpha/beta
singles; `t2_0/t2_1/t2_2` are aa/ab/bb doubles in PySCF occupied/virtual ordering.
The corresponding `l1_*`, `l2_*` are the left-state amplitudes. They are not
independently normalized atomic orbitals. `density_ccsd` is the spin-summed
AO one-body density. `hcore_unperturbed` and `dipole` specify operator
normalization; the exact basis, compact point and state sector are in JSON.
Each NPZ has a SHA-256 link to its source record. No wavefunction files are
loaded with executable pickle serialization.

## Independent comparisons and accuracy limits

First-removal differences are compared with the [NIST evaluated
ionization measurements](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html):
C 11.2602880(11), N 14.53413(4), O 13.618055(7) eV. Displayed rate equivalents
use the common E*/h calibration; eV is only a conversion. The measurements
refer to the lowest physical fine-structure levels, whereas this calculation
omits spin-orbit, scalar relativity, recoil, compact size and radiative terms.
The residual must not be attributed entirely to the solver or called an
uncertainty-level experimental match.

Response is compared separately with the independent **nonrelativistic
CCSD(T)** calculation by [Das and Thakkar, J. Phys. B 31, 2215 (1998)](https://doi.org/10.1088/0953-4075/31/10/011):
approximately 11.67, 7.26 and 5.24 for C, N and O. This is theoretical
cross-validation, not three new experimental measurements. The source's
reported angular components are also retained. No result was tuned to these
values. Different high-level calculations and recommended physical values
need not coincide; these source comparisons do not remove model uncertainty.

A sub-1% last basis change is an initial numerical response/size gate, not
helium-level precision. Energy, response and radial density have separate
convergence histories. Do not infer chemical bond accuracy from accurate
total atomic energies or add isolated-atom correlation energies to an
uncorrelated molecular energy as if that supplied shared correlation.

## Completed campaign and retained limits — 2026-09-16

The final export audits **217 state/perturbation records, 61 zero-field state
exports and six independent small-space FCI controls**. All three campaign
groups completed. The selected C/N states use double-augmented QZ; O uses
triple-augmented QZ. No target value was fitted.

| Quantity | C | N | O |
|---|---:|---:|---:|
| First-removal rate equivalent (PHz) | 2.717794 | 3.512880 | 3.273334 |
| Mean envelope distance, CCSD density (pm) | 62.87038 | 55.68013 | 50.63514 |
| Static angular-term average χ | 11.69706 | 7.26819 | 5.24022 |
| Last cardinal change in average χ | −0.3731% | −0.2565% | −0.6336% |
| Last successive-diffuse change in average χ | +0.08050% | +0.1597% | +0.00594% |
| Largest angular-branch cardinal change, magnitude | 0.4428% | 0.2565% | **1.018%** |
| Last cardinal change in removal interval (eV conversion) | +0.03589 | +0.03291 | +0.09604 |
| Calculated minus NIST removal interval (eV conversion) | −0.02039 | −0.00602 | −0.08063 |

**Grade:** C and N pass the initial static-response convergence checks,
including both angular branches. O's average and successive-diffuse check
pass, but its parallel branch remains cardinal-basis-sensitive. The full O
directional response is **not cleared**. Mean-radius last cardinal changes
are below 0.2% for all three, within this CCSD density definition. Removal
energies still move by 0.03–0.10 eV across the last cardinal step; they are
finite-basis benchmarks, not certified basis-limit or spectroscopic values.
The small-basis FCI removal discrepancies (about 0.002–0.0025 eV) cannot
cancel or bound those larger-basis differences. No combined experimental
uncertainty is claimed.

Next atomic refinement is the O axial cardinal space and the removal-energy
basis/correlation limits; no such next calculation is implied to be running
by this export. The atomic package is usable for explicitly graded controls,
not a declaration that every H/C/N/O input is ready for a long molecular run.

**Molecular handoff:** use these states as benchmarks and carefully mapped
initial representations. Retained-support and additive-energy invariance,
joint configurations, basis-superposition effects and gradients at fixed
geometry must pass in the actual molecular representation before a long
optimization. Static response alone is not a dynamic spectrum or a full bond
curve. The long Molecules run has not been authorized as completed by this file.

## Reproduce

Requires Python, NumPy, SciPy and PySCF 2.14.0. In this directory:

```sh
python correlated_campaign.py --out RESULTS/main --workers 3
python correlated_campaign.py --out RESULTS/diffuse --levels 3 4 --augmentations 2 --workers 2
python correlated_campaign.py --out RESULTS/third-diffuse --symbols O --levels 3 4 --augmentations 3 --workers 1
python correlated_controls.py --out RESULTS/fci
python correlated_atoms_build.py --roots RESULTS/main/records RESULTS/diffuse/records RESULTS/third-diffuse/records --controls RESULTS/fci --out RESULTS/export --strict
```

Run the three campaigns sequentially on memory-limited machines. Saved matching
records are reused. A campaign lock prevents duplicate coordinators for one
output directory. Source hashes, individual logs and state checkpoints permit
an interrupted campaign to be inspected and resumed; absence of a complete
record is never treated as a successful calculation.

Implementation references: [PySCF coupled cluster](https://pyscf.org/user/cc.html),
[SCF and open-shell state treatment](https://pyscf.org/user/scf.html).
