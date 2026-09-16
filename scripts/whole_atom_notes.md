# Whole-atom preparation: correlated helium and response inputs

15 September 2026. Primary exhibit: [The Elements / Atomic Ledger](../poams-periodic-table.html#whole-atom-helium).

## Result in plain language

The accepted measured alpha now calibrates an explicitly specified atomic
energy account. Solving both mobile entries together produces a helium
envelope, binding/removal energies, recoil, angular correlations and
directional response without fitting helium atomic energies or sizes.
The response supplies leading neutral-atom attraction coefficients.

This is a useful effective-model calculation. It is **not** a derivation of
the compact component's internal structure, of the isotope-count assignment,
or of the complete effective operator from the NR axioms. It uses the same
inverse-distance operator already printed in The Molecules, whose mathematics
is the nonrelativistic helium Hamiltonian. There is no hidden new POAMS-only
claim for familiar mathematical machinery.

This calculation was not blind: known helium benchmarks were read before
the computation. They were used for comparison, never as objective-function
targets or coefficients. The preregistered calculation specification is
preserved in the private review record; all executable mathematics is here.

## 1. Inputs and conventions

From the shared 2022 CODATA calibration: alpha = 0.0072973525643, light
inertia m = 9.1093837139e-31 kg, h and c. The newly supplied compact
helium-4 inertia is M = 6.6446573450e-27 kg. Its standard uncertainty is
2.1e-36 kg. The source convention calls this the alpha-particle mass;
this is a measured-inertia correspondence, not a new particle ontology.

Also specified: effective eye strength Z=2, two mobile entries, spatially
symmetric paired-intrinsic singlet sector. The compact I=0 state is a
measured assignment. The state selection is not inferred from the proposed
2-depth : 2-pattern count assignment. Opposite intrinsic slots do not reverse
the mobile/mobile pair sign.

Define K = alpha*hbar*c, a_* = hbar/(m*c*alpha), E_* = m*c^2*alpha^2,
and eta=m/M. These give a_* = 52.91772105 pm and E_*/h = 6.579683920 PHz.
The conventional E_* conversion is 27.21138625 eV. Alpha is not tuned.

**One input to the comparison, not the solve:** the NIST compact rms
charge-distribution radius 1.6785(21) fm. Compared with the calculated
57.81935 pm envelope rms distance, this gives approximately 34,447:1.
Its relative uncertainty from the compact measurement alone is 0.13%.
This is a comparison of differently weighted spatial observables, not a
derived disk-radius or thickness ratio. The point-eye solve never uses this
radius.

## 2. The external whole-atom account

Remove overall translation and use r1,r2 relative to the compact component.
In E_* and a_* units,

    H = -(1+eta)(laplacian_1+laplacian_2)/2
        - eta gradient_1 dot gradient_2
        - 2/r1 - 2/r2 + 1/r12.

This follows from the **specified** pair potential and ordinary three-mass
kinetic energy after removal of the centre of mass. It is not obtained
from counts alone. The positive kinetic form is

    T_mobile = (1/2) integral (|gradient_1 psi|^2+|gradient_2 psi|^2),
    T_recoil = (eta/2) integral |gradient_1 psi+gradient_2 psi|^2.

T_recoil is compact translational recoil, not intrinsic core energy.
The two attractive terms and one mobile/mobile repulsive term are each
included once. The compact internal binding and rest energy are outside
the external energy zero; they do not become zero by being omitted.

The reported finite-mass result, rounded to its useful numerical precision:

| Contribution | E/E_* | E/h, PHz |
|---|---:|---:|
| Mobile kinetic sum | +2.9028848 | +19.10006 |
| Compact translational recoil | +0.000419750 | +0.00276182 |
| Joint compact/mobile attraction | -6.75230634 | -44.42804 |
| Joint mobile/mobile repulsion | +0.945697236 | +6.22239 |
| Total | -2.9033045525 | -19.1028263 |

The virial identity 2T+V=0 is checked independently. There is no physical
necessity to allocate half of every joint term to each region. A chosen
partition of joint energy is a convention, not a measured regional clock.

The one-entry ion has energy -2/(1+eta). Subtracting the neutral energy
gives the **first-removal interval**, 5.94526225 PHz (24.58762905 eV).
This is not the total binding of both mobile entries, about 79.0029416 eV.
NIST gives 24.587389011 eV for the first interval; the residual is
+0.00024004 eV, or +9.76 ppm. No uncertainty-level agreement is claimed.
Relativistic, radiative and finite-size terms are outside this operator.

The independent shared-shape variational control gives only 23.06274 eV
for the first-removal interval. The correlated improvement changes the
joint state, not the pair coefficient or measured inertia.

## 3. How the joint state is computed

Use Hylleraas scalar coordinates s=r1+r2, t=r1-r2 and u=r12. The basis is

    exp(-zeta*s) * (zeta*s)^l * (zeta*t)^(2m) * (zeta*u)^n,
    l+2m+n <= p.

Each basis element is rotationally scalar and exchange-symmetric.
The scale zeta minimizes the calculated energy in a fixed initial bracket
[0.8,3.0]. A nested degree sequence p=2,4,6,8,10,12 is computed; p=14 is
an independent convergence check. Gram conditioning can change the retained
span, so strict nestedness is not asserted after support truncation.

For nonnegative perimetric x,y,z,

    r1=x+z, r2=y+z, u=x+y,
    d^3r1 d^3r2 = 16*pi^2*r1*r2*u dx dy dz.

The exponential factor is exp[-2*zeta*(x+y+2z)]. Product Gauss-Laguerre
integration calculates overlaps, the gradient kinetic form and pair
potentials. The code works in scaled variables; T scales as zeta^2 and
V as zeta. Independent quadrature orders are checked. Radial profiles
use a separate radial/angle integration, not a histogram of those nodes.

The diagonal-normalized Gram eigenvalues define retained support with
relative cutoff 1e-13. At p=12 the scalar basis retains 176 of 252
functions and the directional basis retains 290 of 455. Checks at
1e-12 and 1e-14 disclose sensitivity to this choice. There is no claim
that nearly dependent numerical directions are physical new states.

The infinite-inertia result is -2.9037243718 E_*, within 5.3e-9 E_* of
Schwartz's high-precision same-operator value. Basis/conditioning tests
support about 1e-8 E_* for the energy and 3e-6 absolute for static
response and C6 here. These are convergence indicators, **not** rigorous
error bars or experimental uncertainties. This is double-precision
variational calculation, not the earlier 80-digit alpha substitution.

**Export normalization.** The numerical Gram measure suppresses the common
16*pi^2 angular factor. Therefore the normalized physical spatial amplitude
is zeta^3/(4*pi*a_*^3) times the exported polynomial and exponential, evaluated
at scaled coordinates rho_i=zeta*r_i/a_*. Its full six-dimensional norm is
one. The handoff JSON states this conversion and the directional-response
normalization explicitly; raw polynomial coefficients must not be mistaken
for a separately normalized radial orbital. Spectral strengths already
include the proper normalization and are one-Cartesian squared matrix
elements of D/a_*.

## 4. What “size” and “angular momentum” mean here

The calculated marginal envelope has mean distance 49.19273 pm, rms
distance 57.81935 pm, and mean inter-entry separation 75.26209 pm.
Its one-entry radial distribution integrates to one; the 90% radial
quantile is about 89.81 pm. None is a rigid shell wall or a disk radius.

For the compact-relative angular operators, the scalar state has total
L=0, but separate variances are about 0.0094754*hbar^2 and the cross
correlation is the negative of either variance. Thus

    <(L1+L2)^2> = <L1^2>+<L2^2>+2<L1 dot L2> = 0.

This is a variance/correlation identity, not a pair of classical vectors
observed pointing in opposite directions. The specified intrinsic singlet
and measured I=0 compact state give J=0. Neither that state selection nor
the physical internal core is derived from the isotope-count ratio.
The earlier circular pair's condition L=n*hbar does not apply to this
stationary scalar ground sector.

## 5. Directional response: a necessary molecular input

Apply the weak dimensionless perturbation F*Dz, with Dz=r1z+r2z. A
compensating drive acts on the compact component; the relative deformation
is not a whole-atom thrust prediction. Construct complete L=1 multiplets:

    even-t polynomial * exp(-zeta*s) * (r1+r2),
    odd-t polynomial  * exp(-zeta*s) * (r1-r2).

Both preserve exchange symmetry. The three Cartesian components are
treated equally; their orientation-averaged integral tensors are exact
for a scalar reference state. This angular family is distinct from
product quadrature moments and from a larger occupied scalar basis.

Let A=H_response-E0*S_response and b=<response|Dz|ground>. Then

    chi = 2*b^T A^(-1)*b,
    Delta E/E_* = -chi*F^2/2,
    <Dz>/a_* = -chi*F.

Energy denominators are positive in the retained support. The complete
dipole-strength sum reproduces <Dz^2>; its energy-weighted counterpart
approaches 1+2*eta. An independent finite-field block diagonalization
checks the energy curvature.

Reported chi: infinite eye 1.38319022; finite helium-4 1.38380803.
Independent same-model values are approximately 1.383192174455 and
1.383809986408, respectively. The **physically corrected** helium value
1.38376078 is a different comparison: the present operator omits its
corrections. Do not claim to have calculated them because the numbers
are close.

**Exact H control.** With psi0 proportional to exp(-r), the first-order
correction per unit F is

    delta psi/F = -(r+r^2/2)*cos(theta)*psi0.

It satisfies (H-E0)(delta psi/F) = -z*psi0 and gives chi=4.5 exactly.
Occupied-only space gives zero; one radial r*exp(-r) directional family
gives 4; adding r^2*exp(-r) recovers 4.5. Thus a perfect occupied state
does not establish an adequate molecular response space.

## 6. The long-range coefficient

The response eigenpairs supply positive gaps Delta_n and one-Cartesian
strengths s_n. They represent the resolvent, including continuum response
in a finite pseudostate approximation. They are **not** certified
spectroscopic excitation energies.

    chi(i*w) = 2 sum_n Delta_n*s_n/(Delta_n^2+w^2),
    C6(A,B) = 6 sum_nm s_An*s_Bm/(Delta_An+Delta_Bm).

An independent frequency integral (3/pi)*integral chi_A(i*w)chi_B(i*w) dw
checks the double sum. The imaginary argument is an integration device,
not a physically imaginary clock.

In E_* a_*^6 units:

| Pair | Calculated | Independent benchmark |
|---|---:|---:|
| H/H, infinite eye | 6.49902670539 | 6.49902670541 |
| H/He, infinite eye | 2.82134176 | 2.82134392 |
| He/He, infinite eye | 1.46097620 | 1.46097784 |
| 4He/4He, finite eye | 1.46212122 | mass-corrected output here |

This yields U/E_* approximately -C6/(R/a_*)^6 in the non-overlapping,
nonretarded leading-dipole regime. It is not a complete molecular potential:
short-range repulsion/exchange, higher multipoles and additional corrections
remain. Do not add C6 again if the full correlated molecular solve already
includes the same response.

## 7. Handoff and the next gate

Reusable now: exact H response checks; correlated two-entry energy and
recoil checks; the helium marginal; static/frequency response coefficients;
long-range H/H, H/He and He/He controls. The JSON also exports the retained
state coefficients, basis definitions and response strengths so a later
implementation can be compared rather than trusted by name.

Not supplied: C/N/O many-entry response tensors, a physical planar-core
functional, universal isotope regional ratios, molecular bond lengths,
reaction barriers, or any apparatus propulsion response.

Before a long Molecules search, the corresponding C/N/O angular response
spaces and joint-configuration representation must be recovered and checked.
The earlier water representation/energy-shift/support and fixed-geometry
derivative gates still apply. This helium calculation neither reruns nor
retroactively validates those failed tests. The necessary next step is a
bounded atomic/response and fixed-geometry check, not a multi-day optimizer.

Run `python3 scripts/whole_atom_build.py` with Python and NumPy. It produces
`whole_atom_results.json`, `whole_atom_handoff.json`, the radial CSV and
browser display data. No network access or reference-value fitting occurs
in the solver. The 56 checks are numerical identities, convergence tests and
comparisons—not 56 independent physical discoveries.

## Primary references

- [NIST CODATA table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt): calibration and compact reference radius.
- [NIST helium ionization value](https://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=2): observable comparison.
- [NIST helium state and isotope data](https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable1.htm): I=0 and ground-state correspondence.
- [Schwartz, 2006](https://arxiv.org/abs/math-ph/0605018): independent infinite-inertia energy benchmark.
- [Puchalski, Jentschura and Mohr, 2011](https://www.nist.gov/document/blackbody-radiationpdf): infinite-inertia response benchmark.
- [Puchalski et al., 2020](https://arxiv.org/abs/1912.12242): finite-inertia and physically corrected response comparisons.
- [Yan, Babb, Dalgarno and Drake, 1996](https://arxiv.org/abs/atom-ph/9607002): long-range coefficient benchmarks and response formulation.
