# All-element atomic references — method and reuse

This is an individually solved **mean-field reference campaign for neutral
elements**, not 118 copies of the two-entry helium solution and not a claim
that all isotope cores have been derived. The existing exact hydrogen and
correlated, finite-inertia helium-4 calculations remain higher-grade controls.

## Inputs, model, and units

Measured alpha = 0.0072973525643; the same measured light-entry inertia, h and c
as the Atomic Ledger. Use a* = hbar/(m c alpha), E* = m c^2 alpha^2. In SI,
a* = 52.917721054674 pm and E*/h = 6.579683920428908 PHz.

Each calculation contains Z mobile entries and an infinitely massive point
compact component with strength Z. The selected spin sector is supplied by
the tabulated configuration in PySCF 2.14.0's CONFIGURATION table and a Hund
assignment. This is a state input, not a derivation of the periodic ordering,
isotope spin, or compact-component geometry. Open-shell states can be
anisotropic. The radial profile is the orientation average.

The energy functional is unrestricted Hartree–Fock, with a spin-free
one-electron X2C kinetic/compact-attraction operator and the bare two-entry
Coulomb and exchange terms. The exchange term implements the assumed
antisymmetric state rule; it is not an extra physical interaction. The
operator is conventional-equivalent effective mathematics, consistent with
the existing molecular calculation. It is **not independently derived from
NR ontology by this campaign**.

References: [PySCF SCF documentation](https://pyscf.org/user/scf.html),
[spin-free X2C implementation](https://pyscf.org/_modules/pyscf/x2c/sfx2c1e.html),
[configuration inputs](https://pyscf.org/_modules/pyscf/data/elements.html),
[basis interface](https://pyscf.org/user/gto.html).

## Numerical representation

All-electron Dyall v2z and v3z Gaussian bases, with two additional diffuse
functions in each angular channel through the highest occupied l + 1 (at
least p). The added exponents are the channel minimum divided by 3 and 9.
PySCF 2.14.0 removes overlap eigenvectors with eigenvalues at or below 1e-6
by default. This numerical rank cutoff can leave fewer retained orbitals than
AO functions; the saved coefficients are then rectangular. The verifier
independently checks the overlap spectrum, retained dimension, orbital metric
and entry count. This reduction is disclosed, not a full-basis convergence test.

No observed atomic radius, energy or response is fitted. Basis exponents
are imported numerical representation data, not POAMS-derived constants.

The RI-JK auxiliary basis is generated with an even-tempered ratio of 2.
It factorizes the same pair integrals; it is not a fit to physical data.
Every state is checked against direct four-centre integrals. The resulting
energy and gradient differences are retained in the raw record. Stationarity
is certified for the RI-JK functional actually solved, within the retained
orbital space. No overlap-cutoff convergence claim is made. The direct four-centre
gradient is separately flagged when it exceeds the same 3e-5 threshold; many
references do not pass that direct-integral stationarity test. A small
relative total-energy correction is not proof that the state or its response
is converged with respect to the auxiliary representation.

A spherical radial reference is optimized with the specified angular
configuration; integer spin occupations then seed the unrestricted solve.
The triple-zeta solve starts from the projected double-zeta state. Neither
SCF convergence nor this initialization establishes a global minimum over
all configurations and multiplets. Any change of dominant angular counts
is exported alongside the input counts, not concealed.

## Outputs

* The total **external** binding account and its one-entry, direct and
  exchange contributions. The compact component's internal energy is absent.
* The normalized angular-average radial distribution. Physical large and
  small components are reconstructed using the scalar X2C transformation;
  the small-component density contains the spatial-gradient contribution.
  Mean, rms and 90% radii are different statistics, not shell walls or planar
  core radii. Compact recoil and isotope shifts are not included.
* Occupied orbitals, orbital energies, and state vectors in both spin sectors.
  The negative highest occupied eigenvalue is a **frozen-orbital removal
  estimate**, not a relaxed ionization calculation or measured interval.
  A separate calculation also relaxes the +1 state after removing the highest
  occupied entry from either spin account. Its energy difference is the
  separately labelled Delta-SCF removal interval. The lower converged result
  among these seeds is used; the campaign does not prove the global ionic
  ground multiplet or include correlation corrections.
* Determinant expectations of L^2 and S^2. The imposed value of 2M_S, spin
  contamination and angular-count comparison are disclosed. These are not
  assertions that a determinant is an eigenstate of total J, and they do not
  determine an isotope's intrinsic compact spin.
* Static dipole-response tensors to a small directional perturbation. The
  relaxed response solves the coupled mean-field response equations. Linear
  residuals, symmetry and positivity of every tensor eigenvalue are checked;
  failed responses are withheld. A positive orientation average does not
  rescue a negative directional response. The one-entry case has no pair-response kernel.
* Separately labelled **frozen-mean-field** transition strengths and C6 values.
  Their spectrum holds the self-consistent potential fixed. It must not be
  combined with the relaxed static response as though they were one spectral
  model. These C6 values are reference estimates, not correlated dispersion
  predictions or molecular binding curves.

For a dimensionless perturbation F in direction d, the relaxed tensor gives
Delta E/E* = -F^2 chi_dd/2 and <D_d>/a* = -chi_dd F. The displayed scalar is
Tr(chi)/3, an orientation average of the chosen determinant. No centre-of-mass
thrust follows from an internal displacement.

## Verification and precision

The raw records include entry count, orbital orthonormality, SCF gradient,
radial normalization and nested radial-grid comparison, energy-account
closure, direct-integral check and coupled-response residual. Small-field
energy curvature independently checks the H/He/C/N/O response implementation.
The removal audit separately checks all occupied spin-removal seeds, their
entry/spin bookkeeping, solver-history energies and gradients, and selection
of the lowest stationary tested ion. A positive interval alone is not this
stationarity check or proof of a global ground multiplet.
Small perturbations can expose near-degenerate state rotation and subtraction
error; raw control results are retained rather than averaged away.

The saved larger-basis Y, Zr and Pa determinants, and both basis levels of Lr,
have negative response directions. A separate check rotates the saved orbitals without changing their
occupations: both signs lower the RI-JK energy, confirming instability in a
tested dipole-coupled mode. Their stationary references are retained, but the
relaxed responses are withheld; these are not negative physical ground-state
polarizabilities. The [diagnostic records](all_element_response_diagnostics.json)
include source/state hashes and all three rotation sizes. The coarsest Zr
curvature differs by 10.5% from the linear-response derivative; successively
halved rotations reduce the difference to about 1%, with the energy-lowering
sign unchanged. Total-energy subtraction limits the smallest-step precision.
For Lr the orientation average is positive at both bases, but two tensor
eigenvalues are negative. The larger-basis independent rotation curvature is
-0.00420202 E*, with all tested signed rotations lowering energy and the finest
finite curvature within 0.375%. These are retained-reference instabilities,
not physical ground-state polarizabilities. No symmetry-broken replacement
has been substituted.

Pa/v3z has directional response -849.67956 and normalized rotation curvature
-0.000155072 E*. Its initial smallest-step total-energy subtraction failed
the 2% comparison (15.17% difference), despite both signed energies falling.
Subtracting energies near -27220.9 E* to resolve a change near 5e-10 E* is
cancellation-sensitive. A second evaluation uses the exact quadratic UHF
functional increment, dE = Tr(F0 dD) + Tr(dD V[dD])/2, with SVD sine/cosine
orbital rotations to construct dD without subtracting occupied projectors.
At the same three rotation sizes, curvature differences are 4.45%, 1.23%
and 0.421%; the unchanged finest-step 2% criterion passes. Both signed
increments are negative. The original subtraction results remain in the
record; the response remains withheld. This diagnoses the saved reference,
not a replacement physical ground state. Conversely, a positive
response tensor elsewhere is not a proof of stability against every possible
orbital variation.

The initial response jobs used block Krylov with GMRES residual refinement;
later jobs use symmetry-preserving MINRES with the same GMRES fallback and
the same independent residual acceptance criterion. Saved-state controls
compare the two numerical routes, including an uncertified negative-response
case. No target value or physical term is changed by this solver optimization.

Both basis levels are compared where stationary solutions are available.
Their difference is a sensitivity measure, **not a rigorous error bound**.
Response differences over 5% are visibly flagged as basis-sensitive. The
numerical solver tolerances do not establish experimental precision.

Missing correlation, spin-orbit, two-entry relativistic picture change,
Breit/radiative corrections, finite compact size and recoil remain physical
limitations. High-Z references require these effects for accurate comparison
with experiment. A complete numerical reference across the table does not
mean that all measured atomic observables or core/envelope partitions have
been derived.

## Completed collection and retained exceptions

All 118 elements have stationary RI-JK references at both basis sizes (236
neutral states), and all 118 larger-basis removal intervals are positive for
the separately relaxed tested ions. Of the selected larger-basis states, 114
pass the stated dipole-response checks. Y, Zr, Pa and Lr do not; their relaxed
responses are withheld, not interpolated. Positive dipole-response checks do
not establish stability in every orbital mode or physical ground-state accuracy.

The full saved-state audit passes 7,221 numerical/bookkeeping assertions.
Direct-integral stationarity is a separate limitation: 202 of the 236 states,
including 100 of the 118 selected states, exceed the 3e-5 direct-gradient
threshold. The largest direct gradient is 0.181997 and the largest absolute
energy correction is 0.0172835 E*. These references require auxiliary-integral
refinement before direct-integral precision can be asserted; none was silently
replaced by the direct energy alone.

Ti's relaxed orientation-average response changes +88.94% between bases;
Lu's changes -6.32%. Both are visibly marked basis-sensitive. Mt/v3z changes
its dominant occupied (s,p,d,f) counts from the supplied (14,30,37,28) to
(15,30,36,28), at unchanged specified spin projection. The two mixed
s/d orbitals have leading s weights only 0.500354 and 0.500175, so this
integer count change is a dominant-label threshold crossing, not proof of a
discrete occupation transfer or the physical ground configuration. Db, Ds
and Rg/v3z each retain 384 orbitals from 385 AO functions after the disclosed
overlap cutoff. All three metric/count/projection checks pass, but no
cutoff-convergence claim is made. Larger and smaller bases need not be nested;
energy differences are not guaranteed to be negative.

C/N/O have positive checked tensors, under 0.14% two-basis response differences,
and larger-basis direct gradients below threshold. Their exported states are
usable mean-field starting references and implementation controls, not
correlated precision inputs. Original exact H/correlated He data and methods
are unchanged. No long molecular run is part of this campaign.

## The four retained instabilities and their physical near-degeneracy

Y, Zr, Pa and Lr are the four selected states whose relaxed dipole responses
are withheld from the display. This subsection extracts what the diagnostic
record in
[`all_element_response_diagnostics.json`](all_element_response_diagnostics.json)
says about each, and identifies the physical near-degeneracy each
instability points at. Nothing in this subsection changes a numerical
result; it interprets what is already in the record.

The response-tensor eigenvalues at the retained references (in the calculated
units of the RI-JK response solve, one eigenvalue per Cartesian axis of the
determined-orientation reference):

| State | λ₁ | λ₂ | λ₃ | Neg. dirs |
|---|---:|---:|---:|:---:|
| Y (Z=39) / v3z | −705.28 | +226.18 | +226.18 | 1 |
| Zr (Z=40) / v3z | −4726.62 | −4726.62 | +147.12 | 2 |
| Pa (Z=91) / v3z | −849.68 | +243.38 | +322.92 | 1 |
| Lr (Z=103) / v2z | −18.72 | −18.72 | +118.40 | 2 |
| Lr (Z=103) / v3z | −23.17 | −23.17 | +118.49 | 2 |

At each of the three sampled orbital-rotation norms the finite-difference
energy curvature agrees in sign with the linear-response eigenvalue to
better than the campaign's 2% acceptance threshold, and independently
confirmed for Pa by an exact-quadratic UHF increment (see the paragraph on
Pa/v3z earlier in this file).

What each instability is a signature of:

- **Y (5s² 4d¹, ²D):** A single negative direction. The 4d¹ open shell fixes
  one specific m_l projection out of the fivefold-degenerate d manifold;
  the RI-JK reference thereby breaks the atom's spherical symmetry. The
  negative eigenvalue is a rotational mixing of that occupied 4d with the
  singly-occupied-partner 5s (a well-known 5s/4d near-degeneracy in group
  3). A physical relaxation would state-average across the ²D manifold or
  step up to CASSCF over the (5s,4d) active space.
- **Zr (5s² 4d², ³F):** Two nearly-degenerate negative eigenvalues,
  −4726.62 each. The 4d² occupation fixes two of the five 4d orbitals and
  breaks rotational symmetry in two of the three Cartesian directions.
  This is the classic multireference d-block problem: the ³F ground term
  and low-lying ³P/¹G terms mix strongly, and a single Slater determinant
  cannot represent that mixture. State-averaged multi-reference treatment
  is the standard remedy in the extant literature.
- **Pa ([Rn]5f² 6d¹ 7s², ⁴K₁₁/₂):** One negative direction, moderate
  magnitude (−849.68). Pa carries three open shells (5f, 6d, 7s) and its
  ground multiplet is dense with low-lying alternatives; the near-degenerate
  5f/6d/7s system is textbook-hard even in fully relativistic multi-reference
  calculations. Our scalar-relativistic single-reference solve is already
  outside its guaranteed-adequate regime for actinides of Pa's configuration
  type; the retained instability marks that fact, not an implementation bug.
- **Lr ([Rn]5f¹⁴ 7s² 7p¹ predicted, or 6d¹ 7s² nonrelativistically):** Two
  nearly-degenerate negative eigenvalues at *both* basis levels
  (v2z: −18.72; v3z: −23.17). Lr is a well-studied borderline case where
  the ground configuration depends on whether spin-orbit coupling is included:
  scalar-relativistic (which our RI-JK reference is) predicts 6d¹ 7s², while
  full spin-orbit-inclusive treatments predict 7s² 7p¹ with 7p₁/₂ as the
  singly-occupied spinor. Neither is representable as a single scalar
  non-degenerate reference without spontaneously breaking rotational symmetry
  in some direction. The consistent two-fold degeneracy of the negative
  eigenvalues across two bases is exactly the atomic signature this predicts.

**Common origin.** All four are known multi-reference / near-degenerate
cases in standard quantum chemistry. Our single-reference RI-JK scalar-X2C
treatment cannot resolve them by construction; that is a scope statement
about the campaign, not a claim about POAMS. The retained-reference
instabilities are honest signatures of the underlying physics, made visible
by the campaign's dipole-response acceptance test.

**What is not diagnosed by the response solve.** Whether the *neutral*
state is bound (yes; all 118 have converged, stationary references and
positive first-removal intervals), or whether the exhibit's isotope-count
assignments are affected (no; those are separate bookkeeping quantities).
The response check gates only the propagation of these four elements'
relaxed polarizabilities into downstream molecular work; that gate is
working as designed.

**What would relax the four references.** A state-averaged multi-reference
reference (CASSCF over the physically implicated active space for each
element), a spin-orbit-inclusive one-particle framework (X2C plus mean-field
spin-orbit or a proper Dirac–Coulomb treatment for Pa and Lr), or an
ensemble treatment of the degenerate manifold. None of these has been added
to the campaign; that remains a separately-authorizable multi-reference
extension, not a same-day fix.

**Standing status is unchanged.** The withheld responses stay withheld; the
stationary references stay in the record; no symmetry-broken determinant
has been substituted; no correlated or spin-orbit result is claimed.

## Reusing the files for Molecules

`all_element_results.json` and `all_element_table.csv` provide the readable
summary. `all_element_records.zip` retains individual solver outputs,
comparison levels, transition strengths and methodology. The
`all_element_state_index.json` index lists bounded-size state archives, their
SHA-256 checksums and each contained element/basis record. The selector links
directly to the archive containing its atom. Those archives contain compressed
NumPy records with canonical coefficients, energies,
occupations, explicit basis and molecule metadata. Recreate with PySCF 2.14.0;
load NumPy files with `allow_pickle=False`. For each spin, C has shape (nao, nmo) and C.T S C = I_nmo;
occupations count mobile entries. Do not assume nmo = nao when reloading. These coefficients describe an AO representation, not the
helium six-dimensional correlated amplitude. Do not interchange their
normalizations.

Use these as atomic starting states, basis/response tests and approximate
long-range references. Molecular energies must still be computed in a common
basis with their joint terms and convergence checks. Do not simply add atomic
binding energies or infer a molecular minimum from C6. The correlated H/He
benchmarks remain independent checks on the lower-cost campaign.

## Reproduction

Install Python 3.13, numpy, scipy and `pyscf==2.14.0` in an isolated environment.
Then run:

```sh
python all_element_campaign.py --out ./atomic-results --workers 3
python all_element_build.py --source ./atomic-results --dest .
```

The full published instability diagnosis can be reproduced after Y/v3z,
Zr/v3z, Pa/v3z, Lr/v2z and Lr/v3z are present, using the same environment:

```sh
python all_element_response_diagnostics.py --source ./atomic-results --out ./all_element_response_diagnostics.json --states 39:3 40:3 91:3 103:2 103:3 --stable-energy-differences
```

Use `--states Z:LEVEL` to audit other saved negative-response references
without replacing the default Y/Zr diagnosis. The stable-energy option is
required to reproduce the Pa comparison; it retains total-energy subtractions
as a separate diagnostic. Earlier Y/Zr/Lr records used total-energy subtraction
for their acceptance comparison, as identified by their saved fields.

This diagnostic changes no production state. Its finite rotations test the
specified RI-JK functional, not correlated atomic accuracy. The display attaches
its instability label only when the diagnosed raw-record and saved-state
hashes both match.

The campaign checkpoints each element/basis pair and resumes completed files.
New jobs also save the neutral state before the response solve; a restart can
resume that stage without labelling an unfinished response successful.
An interrupted job is not marked successful. No molecular campaign or patent
modification is performed by these scripts.

For pool resizing, an output-directory marker named
`campaign-pause-PID.json` (PID = that manager process) drains its queued jobs
without interrupting an active atomic/removal calculation. Wait for that
manager and its children to exit before starting a replacement pool. Retire
the marker afterward so it cannot affect a reused process identifier.
