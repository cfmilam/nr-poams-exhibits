# POAMS optical-information link — complete conditional engineering specification

**Edition 1 · 30 October 2026**  
**Design class:** finite, actively guided circulation-information link with a prepared relational angular source, guided transport, angular analysis, energy-conserving relative-phase endpoint, source return and complete attempted-record acquisition.

This specification turns the completed paper construction into one engineering reference: what is supplied, where each operation occurs, what it changes, what remains, and what the final record must contain. It specifies a **conditional finite device**, not a fabricated fiber or a measured optical system. Its geometry is in independently defined comparison coordinates; no meter, hertz, watt, optical wavelength, vendor component or fabrication process is assigned without a material and calibration receipt.

## 1. Scope, status and how to use the specification

The design carries the same original signal and its retained reference through the chain

**source manufacture → local angular preparation → launch → 1,024-edge guided transfer → angular analysis → full cold SWAP → relative-phase endpoint → return of the same source → finite acquired record**.

The engineering claim is closed within the declared material, state, native-operation, clock and readout domains. The result is not a passive-fiber passband, intrinsic optical carrier, generic quantum channel below 1% error, indefinitely reusable streaming device or ready-to-order laboratory bill of materials.

The adopted research authority is the final `paper-completion/REPORT.md` in source collection `fiber-optics-integrated-program-2026-10-29`, together with its conjunctive native, raw-record and final-goal conditions. This document is **new exposition of that version**, not an amendment to its operators, qualified input domain or tests. Where earlier module descriptions give a smaller census, older observer, preliminary error formula or author-stage approval flag, the final assembled values here apply. Mechanism details are inherited only where the final composition retains them. In particular, the optional two-use thermal-stock branch is **not** used.

There are three companion layers:

1. **This technical appendix:** human-readable state, operators, schedules, geometry, accounts, readout, comparison and residual material determinations.
2. **Curated engineering data:** [engineering data guide](data/README.md) contains the exact 2,049 guide phase words, 821 prepared cases, input/output IDs, UNA-65 operands, event/courier schedules, contacts and rational resource accounts, wrapped as derivative references. Large JSON tables use `.json.gz` compression. These are not a newly approved executable or physical measurements.
3. **Source-bound replay:** the separately sealed original executable, its external entry hash and its trusted-runtime/no-concurrent-mutation domain remain authoritative for source-bound record replay. The downloadable engineering references describe the exact finite construction; they do not replace that executable.

All source references below are relative to the source collection. [SOURCE-MAP.md](SOURCE-MAP.md) gives sections and SHA-256 receipts. The source map identifies each scientific definition and its version; the archive contains the engineering references listed in §16.

### 1.1 Provenance classes

| Class | What it means in this design |
|---|---|
| Native physical vocabulary | Circulation, angular momentum, co-signed closure, phase, whole completion, retained exchange and reaction accounts |
| Imported mathematical representation | Finite spin/tensor spaces, projectors, unitary/Born calculus, angular addition, spectral estimates and finite concentration mathematics |
| Explicit constitutive/type premise | Positive native contact catalogue; common-coordinate inverse-length rate matching; finite correlated formed-body geometry; bounded completion and local custody; event-indexed thermal law; UNA-65 atomic capture |
| Derived conditional consequence | Guided packet estimate; same-source intertwining; terminal-plus-work energy; phase probabilities; conserved accounts; error and statistical bounds |
| Unestablished empirical qualification | Actual occurrence of the material types, initial joint-state accuracy, clock and digitizer performance, whole-kit independence, SI scales and useful operation |

Calling a support “LCR,” a word “observed,” or a record “sealed” supplies no physical receipt. The physical meaning is the specific operation and retained factors below.

## 2. Parameter and symbol registry

Unless stated otherwise, the following are fixed design values, not fitted measurements. A dimensionless ratio is not an SI calibration.

| Symbol / name | Meaning and fixed value or domain |
|---|---|
| `k0 = κ` | Positive reference circulation/work tariff; independently defined same-account work unit; SI value unassigned |
| `a0` | Positive comparison-length unit; separately calibrated by geometry comparisons; SI value unassigned |
| Reference time | Independently referenced native completion-count scale; no seconds/post value is assigned |
| `h = Pt + 3 Pq = 2I + L1·L2` | Full nine-state dimer tariff: singlet 0, triplet 1, quintet 3 |
| `J=L1+L2` | Dimer total angular momentum in the imported spin-one-pair representation |
| `s, t_m, q_m` | Complete singlet, three triplet and five quintet states; all remain physical |
| `b=ss`, `a=(tt)J=0` | Two-dimer scalar work words; not single-dimer or angular signal states |
| `u=(a,b)`, `v=(b,a)` | Four-dimer scalar code words, equal tariff `2k0`, total angular momentum zero |
| `R=Q=65,536` | Axial and meridian angular-reference bank sizes; this `Q` is not mechanical resolution |
| `Z=2R+1=131,073`, `X=2Q=131,072` | Primary axial-with-source and meridian triplet counts |
| `N=Z+X=262,145` | Primary angular census; full primary plus copy census `2N=524,290` |
| `b_M=20` | Scalar bits representing `M=-N,…,N`; 80 full dimers |
| `n_rot=n_read=512`, `n_port=128` | Finite symmetric-sweep compilation counts |
| `N_z=6,144`, `V=30,720`, `E=30,719` | Backbone sites, full guide vertices, tree contacts |
| `z_source=2,560`, `z_terminal=3,584` | Axial source/capture stations; separation 1,024 edges |
| `S=12,800`, `T=17,920` | Guide registry vertex IDs corresponding to those axial stations |
| `M=512`, `d=512`, `S_run=1,280` | Envelope cells, translated cells and guided updates; notation distinguished from source `S` |
| `s_ax=2/5`, `θ_ax=asin(2/5)` | Axial matching sine and angle |
| `β=π/2`, `ℓ=10^-6` | Leaf phase and nonzero transverse matching angle |
| `w`, 32-bit phase grid | `θ=2πw/2^32`, integer `w`; 2,049 installed guide words |
| `r_v∈[4/5,4/3]` | Stationary support length/rate parameter fixed before preparation |
| `k_v=κ/r_v∈[3κ/4,5κ/4]` | Guide signal tariff; each owned work dimer has tariff `k_v/2` |
| Guide witness | `k_v/κ=3/4+((37v) mod 65)/128`, except IDs 12,800, 17,915, 17,920 are prequalified at 1 |
| `m_v=1+2deg(v)` | Owned signal/work dimers per guide support: 3, 11 or 13; actual total 153,596 |
| `N_v=h(signal_v)+(1/2)Σowned h(work)` | Conserved support load, `0≤N_v≤21`; canonical loaded sector has every `N_v=1` |
| `L0=√3/2` | Guide comparison reference loop length; not an optical wavelength |
| `B_guide=6,553,600κ` | Positive ordinary guide-contact strain/activation envelope |
| Mechanical `Q_G=2^42`, `Q_BD=2^43` | Finite residual-coordinate resolutions of source-factory families; their `B=100Q^4κ` remains paid |
| Outer held skin `1/256` | Cross-body finite-pose allowance; not independent microscopic displacement at every endpoint |
| Native local service radius 16 | Qualified bounded reciprocal-custody domain in comparison coordinates, not continuous travel speed |
| Endpoint `r_A=1`, `r_B∈[4/5,4/3]` | Two three-owner supports; witness `r_B=8/7`, `k_B=7κ/8` |
| Endpoint `n=0,…,31` | Actual requested tag count in five acquired bits; completed count uses six separate bits |
| `q∈{0,1}`, `reverse∈{0,1}` | Optional `π/2` quadrature and sign reversal; positive reverse tag is `15π/8` |
| `ζ_service=10^-15`, `ζ_bulk=10^-8` | Uniform qualified primitive/service errors, not observed manufacture tolerances |
| `d_in=10^-5` | One initial **joint density half-trace** bound before all preparation, including arbitrary reference |
| `q_halt=10^-4` | One entire-program marked-halt allowance |
| `d_dec=10^-5` | One complete classical-record total-variation decoder allowance |
| 674 exact `C` factors | Independent pure observer factors IDs `0…667` and `684…689`; every other factor/purifier belongs to `A` |
| `P*=1,504,551,385,726,580,037` | Final pre-OFFER native post upper bound |
| Fault denominator | `10,000P*=15,045,513,857,265,800,370,000` |
| Strict OFFER cutoff | `16P*+4=24,072,822,171,625,280,596` |
| Observer suffix | 74,886 posts after `16P*`; cutoff-to-final interval is 74,882 posts |
| Final cursor | `16P*+74,886=24,072,822,171,625,355,478<2^65` |
| Total posts / fresh attempt | `17P*+74,886=25,577,373,557,351,935,515` |
| Scalar output census | 1,811 scalar occurrences, on 7,244 distinct output dimers |
| Prepared cases / repeats | 821 cases × 80,000 intended attempts = 65,680,000 |
| Simultaneous indicators / slack | 17,241 and `η=0.01` |

The source collection defines further *local* indices, endpoint ordinals, contact metadata, cold exponents and program-loop variables. They are fully enumerated in the companion field/contact/ROM data and generator references; none is a free fit parameter.

## 3. Full physical state and supplied material types

### 3.1 Dimer space, scalar words and complements

A physical dimer is the entire nine-state object, not a selected two-state approximation. `h` has sectors `0,1,3` with multiplicities `1,3,5`. All source, work, frame, record and service operations specify what happens on their full complements. A wrong work word, quintet, noncode register, extra excitation, incomplete contact or failed cold preparation remains in the physical inventory and accounts.

On a four-dimer scalar bit, `u=(a,b)` and `v=(b,a)` have equal energy and zero total angular momentum. A code-X exchanges them and is identity elsewhere. Code-controlled X/CNOT and finite carry operations retain their named controls and noncode complements. A scalar ternary snapshot has the isometry

```
u ⊗ 00 → u ⊗ 00
v ⊗ 00 → v ⊗ 01
x ⊗ 00 → x ⊗ 10   for every orthogonal noncode x.
```

The original `x` is retained. Outcome `11` is invalid, not another valid source value. Such a snapshot is an admitted native readout operation with a paid positive implementation; it is not a software classification magically acting on an unmeasured body.

### 3.2 What is supplied

The conditional material catalogue supplies:

- finite nine-state circulation dimers and scalar cold/release organizations;
- positive local scalar/angular contacts and their inverse-address partners;
- formed contacts with specified partner/type/port metadata, reference frames, positive activation, outlet and strain accounts;
- complete rigid/correlated born body templates, not independently movable microscopic points;
- common-coordinate LCR supports with exactly matched owned work rates;
- complete 107,666-dimer controller/courier bodies, 4,096 scalar fields, immutable fault stock and bounded native cadence;
- complete finite local heat packs, untouched passive purifiers and a specified collision opportunity roster;
- the separate independent clock and UNA-65 atomic capture law.

A large positive stiffness or resource count does not derive the physical existence of these types. A specimen either meets its independently specified domain or lies outside this device claim. Unlisted matter dynamics cannot be purchased with a small error-budget term.

### 3.3 One input boundary

Before **any** factory work, the complete nonclock state `A` and reference `R_ref` must satisfy

$$D(\rho_{A,R_{\mathrm{ref}}},\sigma_{A,R_{\mathrm{ref}}})=\frac12\Vert\rho-\sigma\Vert_1\le10^{-5}.$$

Exactly 674 pure independent clock/control factors `C` are tensor-independent as qualified. New supports, heads, records, work banks, settings, scratch, heat packs, apparatus clocks other than the exact `C`, and every purifier belong to `A`. There is no second perfect inlet after the source factory, no per-head closeness coupon and no success-selected requalification.

### 3.4 Cold, positive and heat stock

For a raw initially singlet primary of tariff `t`, a six-dimer pack includes the primary, positive-H transpose copy and compensating pairs. Its inlet bound is `6t`, full capacity `24t`, nominal warm coefficient `3e^-32+5e^-96`; release/presentation/guard costs 2,049 posts and conservatively `2049πζ_service`.

For a charged arrival pair of two dimers at `t=k_v/2`,

$$H_c=t(h_1+h_2)+3t(I-\lvert a\rangle\langle a\rvert),$$

has gap at least `t`, maximum `9t`, and nominal bad weight at most `80e^-64`. Its 16-dimer complete module shares the two actual work outputs with the guide roster and adds 14 retained factors. The added preparation inlet is `9k_v`, capacity `75k_v/2`, and 2,049 posts plus the corresponding service error. All 30,719 guide arrival pairs are actually supplied and charged.

An LCR support's three scalar outputs each use a complete 26-dimer factory; the 84-dimer internal support includes their outputs/remainders plus its positive control/ticket. Internal inlet `147κ`, capacity `592κ`, 6,211 preparation posts and `6149πζ_service` do **not** include the separate member feet, root contacts or owned signal/work energy.

For an eligible local material pair `lo,hi` of fixed gap `Δ>0`, the heat primary is two complete dimers with `H_b=(Δ/2)(h1+h2)`. The scalar ground is `s`, and scalar `a` has energy `Δ`. The full involution exchanges `|hi,s>` and `|lo,a>` and is identity on the entire complement. It commutes with complete energy and angular momentum; in particular `|lo,s,ψ>` is fixed for every allowed angular/reference state `ψ`. Eight additional passive-dual/positive-compensation dimers make a complete ten-dimer pack with inlet `6Δ` and capacity `24Δ`.

Every selected occurrence uses the declared local owner and actual pair gap, not an upper envelope substituted into its Hamiltonian. The final construction uses fresh allocated packs; it does not import the optional two-use savings. Thermal opportunities occur only at specified idle/isolation windows with affected generators OFF, including periodic windows after `2^20` owner nonthermal posts. Thermal posts do not recursively create thermal posts. This is an event-indexed finite environment, not an arbitrary ambient continuous bath.


### 3.5 Formed-contact strain and exact finite shape comparisons

Within each qualified born frame, use the inherited reciprocal-length response c_e=c_e(0)/λ_e and binding tariff

U_e=B[(1/λ_e−1)²−1], with detached U_e=0.

The nonnegative binding convention U_e+B, positive activation/outlet stock and complete invalid-sector rejection tariff retain the energy that a signed potential alone could conceal. Rest-frame variants apply the same declared polynomial to residual coordinates after subtracting the body's actual born frame. They do not freely refit the shape response for a new body.

The original ordinary-contact family has a finite coordinate residual alphabet, fixed partner/port/orientation/occupancy metadata, positive lengths, registered handedness and rigid-block relations. Later current G/B/D/connector and endpoint/guide instances supply their explicitly finer Q grids and born correlated templates, not arbitrary continuum motion. The actual rest vector between every printed endpoint pair determines λ and the finite energy difference. An absolute squared displacement bound, such as 192 for new endpoints, pays capacity; it is not a substitute lower-gap proof.

Every material move updates its actual local master/copy geometry records under the same guarded map, with real activation, release, recoil and thermal/outlet terms. Four noncoplanar anchors and retained orientation determine an actual body pose rather than merely a count of bodies. A fixed common rigid change of chart does not create a new physical orientation resource. Independent endpoint drift, changed residual polynomial, wrong partner or unlisted deformation is not covered by the exact correlated-body qualification.


## 4. Source manufacture and angular preparation

### 4.1 One fixed relational seed

Two axial banks have `R=65,536` triplet dimers each, two meridian banks have `Q=65,536` each, and the source is one further axial triplet. The primary counts are `Z=131,073`, `X=131,072`, `N=262,145`. Reference-description axes are `z,z,x,x`; they do not supply an unaccounted absolute physical axis.

Write the reference vector as `|b_seed>=|Z,Z>_z |X,X>_x=Σ a_(J,M,μ)|J,M,μ>`. The actual admitted target is the total-angular-momentum-zero retained purification

`|Ω_seed>=Σ_(J,M,μ) a_(J,M,μ) [Σ_m (−1)^(J−m)|J,m,μ>_primary |J,−m>_copy /√(2J+1)] |M>_scalar`.

There are 524,290 physical angular dimers plus the same original 80 scalar M constituents. Full primary/copy energy is `2Nk0`, capacity `6Nk0`; M contributes `2b_M k0` and capacity `12b_M k0`. The source is part of the joint state from the start. Its reduced triplet state is `I3/3`; it is not appended afterward to an allegedly still-pure independent bank.

The linear mathematical intertwiner `V U=(U⊗I_copy,M)V` holds for scalar primary operations in the admitted invariant total-J sector. It proves joint-reference comparisons. It is not an extra physical encoder or source-dependent purification.

### 4.2 Actual manufacture and factor routing

The complete source factory uses one G word, 33 B/B† words, 16 prepared-bank phases and the actual ready predicate on 564 newly supplied free work bits plus the original M20. The fixed ordering is G, B, followed by sixteen repetitions of

`good phase → B† → transfer Pz to D → prepared-ready phase → transfer Pz back to B → B`,

then the fixed known phase. The first fifteen phase pairs are ordinary π pairs; the final pair is the following fixed matched pair of the adopted shaping word. The exact word is constructively defined; the 6,449,012,345-bit shaping tape is not claimed to have been materialized or executed.

**The matched pair and its fixed phase.** Set α=asin(1/√409), ℓ=15, θ=(2ℓ+1)α=31α. In the ordered good/bad comparison basis, define

cos ϕ = −cot θ cot(2α),

z = cos α[sin α exp(iϕ) sin θ + cos α cos θ],

exp(iφ) = 1−cos θ/z.

Choose ϕ in [0,π] and the displayed unit-modulus exp(iφ). Here ϕ=1.9564457250617036 rad and φ=0.7747218041259991 rad. The time order is the good-subspace phase ϕ, the actual B†, the prepared-ready phase φ, and B; equivalently S_(B init)(φ) S_good(ϕ). Re(z)=cos θ/2 ensures |exp(iφ)|=1 and makes the residual bad amplitude exactly zero in the ideal α reference problem. If cos θ=0, omit the final pair rather than divide by zero.

For the phase convention S_P(χ)=I+[exp(iχ)−1]P and each ordinary pair implemented in that same order, the fifteen ordinary pairs contribute (−1)^15 relative to the usual positive rotation. Define

u = sin α exp(iϕ) sin θ + cos α cos θ,

a_good = exp(iϕ) sin θ + [exp(iφ)−1] sin α·u,

exp(iγ_finish) = [(−1)^15 a_good]⁻¹.

Thus γ_finish=1.183435424731897 rad modulo 2π. This is the explicit final scalar phase correction for the 33-call shaping word, not a fitted phase. The independently retained G growth corrections are already in G, and each B† is the inverse of the actual entire B payload; none is applied twice. For the compact-band approximate B construction, use these same fixed α-reference phases; the truncation/finish mismatch is part of the stated bounded allowance, not a claim that its actual amplitude is exactly 1/√409.

G owns the same 524,290 target angular dimers. Borrowed B and D primaries are subtracted once, not recreated. The original targets are Pz: source `sa.data(1)` then `sa` and `ra` axial frames; Px: `sm` and `rm` meridian frames; copies: the prescribed `sa`/`sm` slots and root factor in their retained coupling order; M: the original root slots grouped four per bit.

A 1,071-vertex, 1,070-edge connector tree joins G/B/D and the four optical annexes without identifying distinct roots. It routes 393,218 initial G→B/D factor movements, 4,194,336 Pz interchanges and 524,370 final placements: **5,111,924 complete actual dimer movements**. Every movement releases the real mounts, carries that complete factor/reaction/reference through named local custody, vacates the source and attaches the destination. It is not a free permutation of semantic labels. Both loaded delivery and empty-controller return are charged.

All transported angular targets have the same born endpoint vector `(1/32,−1/32,1/32)` with squared length `3/1024`. Their mount lengths and body IDs are explicit. The connector's angular neighborhoods and actual free scalar work remain separate. Exact scalar source sequences, source-factory formulas, full body generators and adopted section references are in the source map; the giant finite word is specified by its generator rather than an impractical printed tape.


### 4.2A Constructive source words and coefficient specification

The source program is finite even where it is too large to print as a gate tape. Its constructive definitions are part of the engineering specification, not a request for a later unspecified preparation mechanism.

**G, maximal-bank singlet growth.** For each bank length L, grow F=1,…,L−1, giving 262,143 total steps. The local angular sum C_F=2F+J_new·F_old has comparison eigenvalues 3F,2F−1,F−1. With p_F=(2F+3)/[3(2F+1)], use β_F=acos(1−1/(2p_F)). The exact local fresh-pair phase and grown-top phase map the old-bank singlet plus fresh microscopic singlet to the grown-bank singlet with a known retained phase correction.

On the full two-dimer local contact, P1=L(3I−L)/2 and P3=L(L−I)/6. A nine-dimer mobile body—new angular dimer plus two four-dimer term bits—uses root probabilities (1/2,1/8,3/8,0) and the actual first-F-leaf position loader. Its full leaf operation selects the P1 or P3 reflection at its locally occupied leaf, is identity on absent/noncode/wrong-species sectors, and acts on at most 15 external dimers (21 including its switched mode/ticket). No second coherent address register or remote leaf oracle is supplied.

For the full involution U_F=T†P_F†D_FP_FT, the local root/term00 block is A_F=C_F/(4F). Sixteen actual U_F queries alternate with a 17-phase palindrome. The supplied table [growth-phase-words.u64le.gz](data/growth-phase-words.u64le.gz) contains every row F=1,…,131,072. Each row has twelve little-endian unsigned 64-bit entries: nine symmetric θ words, β, negative γ and middle-wrap parity. Angles use the 60-bit full-turn grid, and θ/β are even. The fixed correction is

γ_word = [(2Σ_(j<8)θ_j + θ_8 + β + 2^59)/2 + parity·2^59] mod 2^60.

The middle-wrap parity is part of the actual finite word. It cannot be dropped as an “irrelevant” global sign when that word participates coherently in the seed. Both bank growths reuse the appropriate actual F rows with their named bank control. Every erroneous/off-root state remains in the full physical operation; the comparison telescopes along ideal prefixes rather than reprojecting a noisy bank.

**B, coefficient shaping.** For the retained phase convention Z=X+1, N=Z+X and 1≤M≤J≤N, the nonnegative coefficient probability is exactly

w_JM = (2J+1)(J−1)!(J+M)!(2Z)!(2X)! /
[(N+J+1)!(N−J)!(M−1)!²(J+1)!(J−M)!4^X],

with recurrence w_(J,M+1)/w_JM=(J+M+1)(J−M)/M². All outside coefficients vanish. Thus the finite scalar compiler uses exact factorial/rational recurrences, not an arbitrary fitted amplitude table or an absolute-value replacement of signs.

The fixed compact bands are J=181269…189460 (8,192 labels), M=129025…133120 (4,096 labels). The target tail union is below 2.614×10^-14 and is charged, not removed by measured postselection. For each J, let p^T_J be the sum over retained M. The conditional M amplitudes are √(w_JM/p^T_J); the actual success-flag rotation has amplitude √[D p^T_J/((2J+1)409)], where D=(2Z+1)(2X+1) is the dimension of the product of the two maximal-spin primary banks. Its success arm is √(τ/409)Ω_T and its reject arm remains physical. The fixed 33-call finish uses α0=asin(1/√409); normalization, discarded-target comparison and finish mismatch have the stated bounded allowance, not a claimed unconditioned exact Ω output.

**Exact clock amplitudes.** Let c_t be the coefficient of x^t in (1+x+⋯+x^(L_clock−1))^r. Its complete finite definition is

c_t = Σ_(k=0)^min(r,⌊t/L_clock⌋) (−1)^k binom(r,k) binom(t−kL_clock+r−1,r−1)

for 0≤t≤r(L_clock−1), and c_t=0 outside that support. With Z_clock=Σ_(t=0)^(T_clock−1)c_t², prepare Σ_t c_t/√Z_clock |t⟩ in the enclosing 25-bit word; all t≥T_clock have zero ideal amplitude. The actual binary preparation tree uses exact integer weights c_t²: at each node, the right-child probability is its weight sum divided by the total node weight (zero if the whole node is zero). `tree_levels`, `theta_interval` and `gray_interval_word` implement its rationally isolated finite rotations, and the full physical operation retains its noncode/unused complement. These are amplitudes of a declared clock-control state, not random samples of an asserted distribution.

The selected constructive shaping clock is **mixed radix** κ_clock=3, b=23, T_clock=25,165,824, convolution rectangles r=12, L_clock=2,097,152, scalar angle precision p=48. This κ_clock is not the work scale κ. Register t=a+2^b h uses b low bits and two high bits. Controlled low powers plus the actual [h≥s] predicates execute T_clock−1 real walks per forward clock block. The finite Fourier factorization is F_3 on h, explicit twiddle phases a·c/T_clock, then F_(2^b) on a; unused binary states retain their prescribed complement.

The F_3 operation is explicitly the following **time-ordered** six-operation word on high-register labels 0,1,2, leaving label 3 untouched:

phase₂(π) → R₁₂(π/4) → phase₁(π) → phase₂(3π/2) → R₀₁(arccos(1/√3)) → R₁₂(π/4).

Here phase_j(χ) multiplies label j by exp(iχ), and R_ab(χ) has ordered a,b block [[cos χ,−sin χ],[sin χ,cos χ]]. Their product is the positive-sign F_3 matrix (F_3)_(j,k)=exp(2πijk/3)/√3, with identity on label 3. The nontrivial angle is isolated by the same rational compiler using probability 2/3; it is not an unspecified Fourier oracle. For output k=c+3d, the exact identity (a+2^b h)(c+3d)/(3·2^b)=hc/3+ac/T_clock+ad/2^b+hd fixes the twiddle sign and digit order. The finite phase tapes, masks and complement extensions remain paid.

Across the 33 B/B† words there are 1,660,944,318 actual walk calls, 1,650 copied physical power blocks and 570,812,820,996 scalar payload pulses, each with its applicable positive switching account. The packet pose precision is 52 bits; it cannot be replaced by the 48-bit scalar angle precision. The complete ideal synthesis upper is 0.0008987847, included in the final whole-source error rather than added a second time.

The immutable two-positive-tape program has exactly 6,449,012,345 scalar bits. It is specified by the finite coefficient/threshold/Gray-angle compiler; it is not claimed to be materialized here. The 584 logical work bits include the original M20 and 564 physically separate supplied work bits. The 16 prepared-ready reflections use all 584 actual inputs, giving 18,688 uses with 78,938,112 paid ready-route posts. The 400-bit outer-phase cartridge has actual 50-bit forward/inverse phase tapes; 3,300 delivered bits use 13,200 full scalar SWAPs. Native history remains forward on every payload adjoint.

The companion readable normative source functions include the actual selected compiler and growth indexing, while the source map locates the complete current derivations and certificate. Source program counts are not extra stock to add twice to the final bill: the final integrated census already includes their physical supply, routing, heat, controls and retained copies.


### 4.3 Local angular operations

For source or receiver dimer `s` and one frame dimer `f`, use `P=Pt(s)Pt(f)` and `A=J_s·J_f`. The two complete local generators are

`L=P(A+2I)`, `Q₂=P(2A²+A−I)`.

They are positive on all 81 two-dimer states, rotational scalars and individually commute with `h_s,h_f`. Their triplet×triplet eigenvalues are respectively `(0,1,3)` and `(5,0,2)`; they are zero outside that sector. The quadratic contact is an explicit constitutive operation, not a claim that bilinear exchange implements its own square.

For a K-dimer co-sense frame, the joint-reference comparison errors are `|t|/√K` for L and `√5|t|/√K` for Q₂. Each finite average is implemented by positive symmetric sweeps, not a free collective pulse. With `n` sweeps their full-space compilation bounds are `(3|t|)^3/(3n²)` and `(5|t|)^3/(3n²)`.

Source commands arrive only after common readiness; they remain local. All choices share a padded horizon and launch address. No command-dependent success, retry or early timing message reaches the receiver before the data contact.

| Source command | Q₂ meridian angle | Subsequent L axial time |
|---|---:|---:|
| +Z | 0 | 0 |
| −Z | π | 0 |
| +X | π/2 | π/4 |
| −X | π/2 | 3π/4 |
| +Y | π/2 | π/2 |
| −Y | π/2 | 0 |
| blank | 0 | 0, with actual full source↔local-singlet dump SWAP |

The blank stores source energy and correlations in its source-side dump. It does not choose a new globally prepared blank state.

Receiver analysis is Z: no rotation; Y: Q₂ angle π/2; X: L time π/4 then Q₂ angle π/2. It then runs the circular instrument, unless READ is deliberately omitted with the same reserved horizon. Six bright inputs plus blank against X/Y/Z give 21 ordinary angular cases; three +Z/omit-READ cases complete 24 original controls.

### 4.4 Circular instrument and angular retention

The receiver uses its finite axial bank, the scalar phase circuit on two four-dimer bits, classification, then the **actual inverse** phase circuit. The collective comparison spectrum is `3, 2−1/R, 1−1/R` on the qualified maximal-frame branches, zero on singlet/quintet data. The circuit uses H on both bits, controlled times π/2 and π, then the positive-sign four-point Fourier transform. It assigns integer phase codes 3, 2, 1; classifier code 3 is plus, 1 or 2 minus, 0 reject.

The full record uses `r0=(a,b,b)`, `r1=(b,a,b)`, `rr=(b,b,a)` on six dimers, and a four-dimer completion word W. The root classifier is a bounded full-space scalar involution on 18 dimers, or 22 including the actual READ bit. Phase inverses use positive `3I−L` and the paid scalar phase corrections, not modulo-2π inversion of a noninteger spectrum. All workspace and reference recoil remain.

The pointer comparison bound is

`√[2/((R+1)(2R+1))]+3π/R+81π³/(4n_read²)`.

Completed plus, minus and reject are `(v,r0)`, `(v,r1)`, `(v,rr)`; omit-READ is `(u,r0)`. The third triplet component is retained and uses the real angular effect, not a fabricated two-level probability. All other records remain unavailable or unexpected.

## 5. Guide geometry and conserved local operation

### 5.1 Finite comb, supports and contacts

There are 6,144 backbone vertices `(0,0,z)` and four leaves `(±1,0,z),(0,±1,z)` at each station. Leaves have no axial or mutual edges. The result is a tree with 30,720 signal vertices and 30,719 contacts, not an isotropic glass cylinder. All four transverse contacts are nonzero. Positional update bands are not polarizations: all three triplet angular components are transported identically.

The tree is rooted at axial station 2,560. At each edge `u→v`, four work dimers form upstream-owned pair A at tariff `k_u/2` per dimer and downstream-owned pair B at `k_v/2`. **Both complete pairs are mounted on the arrival side.** Every support owns its signal plus two work dimers per incident edge, with at most 13 owned dimers. Total owned membership is `V+4E=153,596`.

Each vertex has an 84-dimer scalar support and the paid member endpoints/feet, four root contacts, WORK frame and distinct RIBBON caddy. Frame braces and cross-frame comparisons are counted separately; common names do not identify physical roots. The born body atlas, complete caddy endpoints, head interiors, endpoint ports and finite-pose exclusions are specified by rational generators. One full courier has 107,666 dimers, not merely a plotted anchor. The finite-pose proof does not claim a continuous swept volume through intervening material.

### 5.2 Exact local rate matching

The LCR-m support fixes one common stationary coordinate `r_v`; signal length is `r_vL0`, owned work length `2r_vL0`, tariffs `κ/r_v` and `κ/(2r_v)`. This is an exact material/type relation. It is not obtained by making a spring arbitrarily stiff.

The complete energy is

`H_free=Σ_v[f_v(q_v)+κN_v/r_v(q_v)]+all independent service bodies`.

`N_v=h(signal_v)+½Σowned h(work)` commutes with every permitted edge K, Pt tag, guard and presence operation. On the canonical loaded packet every `N_v=1`; failed/wrong-work/full loading sectors retain their actual `0…21` load. No projection replaces N by 1 in a failed state.

For this guide's ordinary positive envelope `B=6,553,600κ`, maximum angular load is `26.25κ`. A fixed-N residual-shape gap is at least `B−10.5κ`; a cross-sector envelope instead uses `B−26.25κ`. Both exceed B/2 in the admitted finite residual model. The lower-gap proof comes from that model and its fixed rest-frame polynomial, not from body size or a loose upper displacement bound.

If independent work-rate errors produce `d_lock=k_u+g_B−k_v−g_A`, then `||[K,H_actual]||=|d_lock|/2`. Every nonzero `d_lock` is outside exact matching; exact energy pinching would remove the desired transition. There is no open tolerance to arbitrary within-spine drift. Different **stationary spines** in the stated interval are allowed.

### 5.3 Full six-dimer exchange

Define three active columns in each direction:

```
L_m = |t_m,s; b_A,a_B>
R_m = |s,t_m; a_A,b_B>
T_e = Σ_m |R_m><L_m|
K_e = (P_L+P_R−T_e−T_e†)/2.
```

K is a rank-three positive projector on the full `9^6=531,441` physical space. Its kernel—including wrong work, extra loading and quintets—is retained. It commutes with full local energy and all total-AM components in the LCR domain.

The actual phase-flat guided gate is

`Ĝ_e(θ)=exp(+iθPt(u)) exp(+iθPt(v)) exp(−2iθK_e)`.

The two Pt tags are indispensable: the K exponential's local block phase is not a common global phase. Each positive angle is compiled modulo 2π. For a base 32-bit word w the three words are `2w,−w,−w mod 2^32`; this correlated arithmetic gives the intended edge angle. On mixed enable sectors the actual phases remain, e.g. one enabled endpoint can retain its Pt phase even while exchange is disabled.

### 5.4 Retained work is the coherence construction

For vertex x, let `w(x)` be spent exactly on the unique root→x path and ready elsewhere. The isometry

`J|x,m>=|triplet m at x; singlets elsewhere>⊗|w(x)>`

intertwines each local physical gate with `exp(iθX_e)⊗I_triplet`. It holds for any finite sequence and arbitrary retained reference. Reverse crossings restore the traversed edge word; every walk ending at x has the same w(x). Recombination therefore requires no work erasure, projection, new battery or discarded path record.

The **joint** signal/work packet is coherent. Tracing out work generally destroys bare position coherence; no contrary reduced-state claim is made. A successful catcher leaves the common spent word at the terminal. Cycles are excluded: two different paths around a cycle can retain orthogonal work histories at the same position and invalidate the one-word-per-vertex reduction.

### 5.5 Energy and angular momentum at a cut

At a completed OFF crossing,

$$\left(\Delta E_{\mathrm{signal},u},\Delta E_{\mathrm{signal},v},\Delta E_A,\Delta E_B\right)=\left(-k_u,+k_v,+k_u,-k_v\right).$$

The arrival side containing both work pairs receives `k_u`. For a complete rooted subtree its gained energy operator is `k_u Π_subtree⊗I_triplet`; scalar work carries no hidden angular sink. The carried triplet carries its true angular representation, and the source partner remains retained outside.

Support load and good scalar state return at a completed boundary and hence contribute zero **change**, not zero stock. During an active switch the positive contact, F, tickets, interaction, recoil and outlets must also be counted. The simple OFF-cut formula is not an unfinished-switch account.

## 6. Guide loader, update word, catcher and controls

### 6.1 Finite compact packet

Use 512 cells with normalized envelope

$$f_j=\frac{\sin^3(\pi j/513)}{\sqrt{5\cdot513/16}},\qquad j=1,\ldots,512,$$

and spinor `v0=(√0.7,i√0.3)`, packet `ψ_j=(−1)^(j−1)f_j v0`. These two entries are positional components of one full triplet, not two polarizations.

The 1,024 positive site weights alternate `√0.7 f_j` and `√0.3 f_j`. Starting with the actual source at station 2,560, execute 1,023 adjacent gates with

$$\theta_j=\operatorname{atan2}\!\left(\sqrt{\sum_{l>j}w_l^2},w_j\right).$$

The sequential i-phases produce the exact packet. The catcher uses the translated edges in reverse with positive inverse angle `2π−θ_j`. Loaded support is 2,560…3,583; target support is 3,584…4,607. The complete finite backbone has sufficient 1,280-cell guards for the finite support cone; no tail is discarded to obtain a nominal answer.

### 6.2 One update and count

Run 1,280 updates, each in this order:

1. Backbone even matching at `asin(2/5)`.
2. Backbone odd matching at the same angle.
3. All-leaf Pt phase `exp(−iπPt/2)`.
4. Four successive transverse matchings at `ℓ=10^-6`.

Every physical gate executes on blanks as well as bright inputs. A deterministic scheduler never reads the signal's position. The fixed 2,049×32-bit guide ROM contains all loader words and update constants; there are 2,049 separate cold word bodies rather than a renamed enormous head. [guide-rom.json](data/guide-rom.json) is the exact ordered integer array.

### 6.3 Derived guided and transport bounds

Let `b=4sin(asin(0.4)/2)`, `r=(b+√2)/2`, `g=(√2−b)/2`, `a=4ℓ`. For `a<g`, the finite-open guided projector satisfies

`e_P≤ra/[g(g−a)]`, `rank P=6,144`.

At the witness `e_P<5.00924×10^-5`, uniformly for every integer update; bare-core cladding probability is below `4e_P²<1.00370×10^-8`. This is a guide-subsystem comparison, not the final apparatus failure probability.

With `t=e_P²/(1+√(1−e_P²))`,

`δ_eff=2t+a²/2+2ae_P+e_P²`,

and the core-to-comb evolution discrepancy is at most `2√2e_P+1280δ_eff<0.000148630`. The compact packet axial bound is

`ε_ax=3π/[√10(M+1)]+d(1+2s_ax²)π⁶√(369/5)/[192(1−s_ax²)^(3/2)(M+1)³]<0.006089405`.

After the actual 32-bit words and bulk/switch/presence services, the isolated transport core norm bound is `0.00690570077622`. It is not substituted for the full assembled `ε≈0.09618`. The numerical finite-word capture probability `0.9999929107044572` is a synthetic diagnostic, not a measured delivery efficiency or the guaranteed complete-device probability.

### 6.4 Delivered controls and physical presence

- **PASS:** complete loader/run/catcher, then terminal presence, angular port and instrument.
- **Source BLANK:** actual upstream source/dump SWAP. K annihilates a blank even with arbitrary charged work; presence is ideally 0.
- **Guide CUT:** actual delivered CUT word zeroes every occurrence of the compiled axial K at station 2,560, including loader/run/inverse. Positive wrappers/padding and Pt tags remain. No leaf path bypasses the cut. Presence is ideally 0.
- **Guide ROUNDTRIP:** execute the complete finite word followed by its positive inverse in exact reverse order **before** presence. The signal and work return to the source; native histories do not reverse or reset. It is not also delivered to the endpoint.

The terminal presence pointer has four actual dimers and generator `Pt(T)⊗D`, where `D=(|u>−|v>)(<u|−<v|)/2`. A paid positive π operation writes occupancy without reading or replacing angular information. Source, port, angular instrument and all failed complements remain part of one account.

## 7. Same-source endpoint, driven phase and complete return

### 7.1 The exact borrowed source

After the guide port and angular readout, a full equal-tariff **9×9 dimer SWAP** transfers the actual receiver data into the qualified cold receiver `ra,r1,slot24`. On the nominal singlet recipient,

`ρ_BX⊗|s><s|_C → |s><s|_B⊗ρ_CX`.

The post-SWAP C is the endpoint source. The newly emptied B is not. C's immutable ID is `('ra',('d',6,0,2,0))`; its born endpoints are `(89,−77,229441)/64` and `(91,−79,229443)/64`, and its tariff remains `k_A=k0`. Arbitrary source/reference correlations follow this same full factor.

Outgoing transport releases nine real original mounts and installs nineteen at the endpoint. Return reverses those mount obligations; 112 thermal release/return windows are paid. The R1 route has 43 vertices, 42 edges and 85 complete caddies. Source movement uses 24 route edges plus two local pickup/docking guards per direction, charged as 26 units each: **52 priced source units**. No retired old endpoint route is reused by name.

### 7.2 Local converter energy

Spine A owns `S,WA0,WA1`; spine B owns `T,WB0,WB1`. Their stationary rates are `k_A=k0`, `k_B=k0/r_B`, with `r_B=8/7` at the witness. The same full rank-three exchange K acts on scalar-work-dressed L/R triplet columns. Both complete work pairs lie beyond the local arrival cut.

R changes energies by `(-k_A,+k_B,+k_A,−k_B)` in `(S,T,WA,WB)`. At the witness, terminal `7k0/8` plus net work `k0/8` equals the original `k0`. Terminal energy alone is not the delivered-energy account.

A local post-guide blank is a paid controlled full `S↔D` SWAP. It puts the actually arriving source and its reference in D while S/T are singlets; it is distinct from absence of any arriving excitation.

### 7.3 Phase protocol

For FRINGE execute `U(π/2)`, n actual acknowledged tags `exp(−iγPt(T))`, optional quadrature `qπ/2`, then `U(3π/2)`, where `U(θ)=exp(−iθK)`, forward `γ=π/8`, reverse positive angle `15π/8`. The right probability is

$$p_R=\frac{1-\cos\!\left((-1)^{\mathrm{reverse}}\frac{n\pi}{8}+\frac{q\pi}{2}\right)}2.$$

The last operation **recombines/interferes** L and R; it does not erase their coherence or reset work. Free OFF dwell does not create a new relative phase because the complete dressed-column energies are equal. This is driven relative phase of the same signal/work object, not an intrinsic carrier.

| Endpoint setting | Ideal result for an arriving bright signal |
|---|---|
| FRINGE | L probability `1−p_R`, R probability `p_R` |
| CUT | L with probability 1; both exchange actions absent, installed supplies/tags retained |
| NO_FIRST | L/R each 1/2, independent of tag setting |
| TRANSFER | R with probability 1; one positive π exchange, `n=q=reverse=0` |
| Local blank | D with probability 1 |
| No signal from guide/source | EMPTY with probability 1; D remains unoccupied |

Tags alias modulo 16. Reversal leaves the q=0 fringe even and reverses the quadrature sine. These are prospective consequences of the specified finite controls.

For the qualified canonical bright input, endpoint effect pulls back to `p_z I_triplet`. Therefore the ideal joint law with previous angular outcome y is `Pr(y,z)=q_y p_z` even with arbitrary retained angular/reference entanglement. The noisy output need not factorize; its deviation uses the **linear** record-TV bound. Generic third-component inputs use the actual angular effect, not the registered two-bright table.

### 7.4 Endpoint program and return

The finite controller holds the actual blank, quadrature, direction, five n bits, three mode bits, six completed-count bits and two STOP bits. It checks type/arrival, executes blank and selected first mixer, then 31 fixed tag slots. Each tag slot computes `i<n` reversibly, uses the literal enable word, performs or pads the tag, and increments the separate count **only after its completion ACK**. Optional quadrature and final mixer follow. Ten actual flag opcodes complete before STOP becomes `(1,0)`. No requested n or cursor is substituted for the completed count.

The maximum schedule is 45 data invocations: one blank, two exchange, 31 tags, one quadrature, ten flags. Each invocation has a 19-ticket positive bank, 40 positive π switch slots, 32 data posts and 64 guard posts: 1,376 posts. The logical-angle maximum is `573π/8`. The separate padded control program has 4,096 scalar slots (the explicit comparator/carry requirement is at most 3,170), each paid for 32 native posts and its service error. Cold preparation and source movement are additional.

Return releases nineteen endpoint source mounts, restores the original nine mounts and retains the same source, reactions and correlations. No cold reset is required or claimed. The 40 flag constituents are carried by two complete heads in loads 24 and 16. Endpoint head feet are released before motion. Four measurement heads carrying 24/22/24/16 external constituents remain held in the observer bays. They are distinct from the four empty controllers that must complete their home returns. Spent work, source, flags, pointer phase workspace, controllers and heat are retained through the final seal.

## 8. Native controls, causal receipts and immutable custody

### 8.1 Actual controller fields

The guide's complete 4,096-field heads retain disjoint algorithm/native/cargo/fault roles. The working guide program uses 180-bit time/deadline fields, 32-bit phase words and fixed reversible comparison/carry logic; no 128-bit counter is silently substituted for the longer schedule. A 4,096-operation scalar block bounds two 180-bit ripple passes, four 32-bit phase passes, four 128-bit cursor/fault comparisons and the actual route/CUT compilation/audit work. Each scalar slot costs 32 native posts, not one free classical instruction.

All 24 cargo descriptors are 64-bit actual fields. Guide descriptors 0…15 use fields 1815…2838; descriptors 16…23 use 512…1023. These are disjoint from algorithm fields 0…511, reserved 1024…1046 and lifecycle fields. The actual physical lookup appears in the guide-native companion, not in a diagram's label.

The endpoint controller's key fields are:

| Fields | Physical role |
|---|---|
| 0…127 | Native cursor |
| 128…255 | Working fault word; immutable initial copy is separately retained |
| 256…639 | Program |
| 640…767 | Position |
| 768…2303 | 24×64 grip descriptors |
| 2304…2559 | Mailbox |
| 2560…2687 | Inlet |
| 2688…2943 | Audit |
| 2944…3007 | Status |
| 3064…3191 | Endpoint 128-bit actual bay token |
| 3200…3210 | Blank, quadrature, direction, n0…n4, mode0…mode2 |
| 3211…3216 | Actual completed-tag count |
| 3217…3218 | STOP0,STOP1 |
| 3219 | Actual enabled-tag word exposed to its formed opcode |
| 3220…3263 | Iteration, borrow/carry and reversible mode/completion scratch |
| 3264…3391 | Frozen global deadline |
| 3520…3530 | Immutable incoming setting copy |
| Later assigned ranges | Final successor event/scratch/baton/handshake allocations; see acquisition-fields data |

The final successor's field file is controlling where it explicitly allocates previously reserved fields. Guide pre-existing completed/prepared fields 3851…3874 are not overwritten. New final handshake words 3880…3882 remain separate from the freeze of acquisition-visible words. A field declared reserved in an earlier component is not automatically free; final lifetimes and actual constituents govern.

### 8.2 Four physical 128-bit tokens

Tokens use least-significant-bit-first encoding with namespace 0x554e4131:

| Bits | Interpretation |
|---|---|
| 0…31 | Fixed namespace |
| 32…63 | Body number 1,2,3,4 for the actual held head |
| 64…95 | Native head generation 0 or 1; distinct from attempt-envelope generation |
| 96…100 | Actual ACK-counted occupied external constituents |
| 101 | STOP/release completed |
| 102 | Actual occupied reciprocal hold at the assigned bay |
| 103…127 | All 25 bits fixed zero |

Each starts from its original paid cold factory, tests valid zero code and performs paid header X operations. Each literal constituent pickup follows the existing grip/reciprocal-port completion, then increments the token count. Complete cargo lists have sizes 24/22/24/16. STOP/release and occupied bits require the physical completion events, not requested commands. The source token freezes after docking; its final copied bits alone cannot establish actual head identity without the typed custody/port relation.

There are exactly 512 once-only descriptor destination copies **after accepted OFFER**. Additionally, 512 ternary source diagnoses retain the counterexample in which a noncode source leaves a destination apparently zero. Source complement 10, output N and pair 11 remain failure. Destination, original source diagnosis and the full fixed token convention must all agree.

### 8.3 Sixty-four event observands

| Events | Actual completed operation or condition |
|---|---|
| 0…30 | Enabled endpoint tag completion ACK; expected i<actual n |
| 31…40 | Ten flag-operation completions, independent of whether the measured flag is 0 or 1 |
| 41 | STOP operation after all required flag completions/count relation |
| 42 | Same original source restored, acknowledged by actual local receipt path |
| 43 | Endpoint's 16 exposed feet released |
| 44 | Guide's 72 exposed eligible feet released |
| 45 | Guide angular/full cold-SWAP and OFF dependency before endpoint borrow |
| 46…49 | Four distinct empty controllers returned to their actual home berths |
| 50…53 | Four actual occupied bay holds, with valid tokens/cargo |
| 54…57 | Acquisition-visible fields frozen after late return receipts and scratch uncompute |
| 58 | Local preoffer conjunction of real held bays and required completion/custody relations |
| 59 | Actual OFFER accepted |
| 60 | Acquisition started |
| 61 | Acquisition completed successfully |
| 62 | Independent final clock seal |
| 63 | Reserved, never written; expected zero |

Every ordinary event begins zero on the nominal cold branch, has one named completed-address writer and remains immutable afterward. Its full source/output word and physical dependency graph appear in [acquisition-events.json.gz](data/acquisition-events.json.gz). Cold/fault/input errors are retained in the single error account, not eliminated by a logical induction.

The good-path implication comes from the acyclic graph: a receipt is written only at the completed address or after its actual adjacent predecessor packet. Induction establishes that a valid event58 implies the named operations. The converse requires the explicit bounded-progress domain; eventual fairness is not enough to meet a finite cutoff.

### 8.4 Receipt couriers and local packets

The four empty return controllers are guide dga(17921,1), dga(17922,1) and endpoint rsj(1,0), rsj(1,1). They are not the four spent measurement heads. Source return and endpoint return ACKs follow actual stationary neighboring owners. Guide home returns have no intermediate stationary mailboxes: two separately supplied full 107,666-dimer receipt couriers carry those receipts.

Each guide receipt courier makes a 64-hop tour from its staging berth through the original empty-delivery corridor to the exterior port and back. Total: 128 new whole-body hops. The original returned controllers remain at home; couriers never enter an occupied home berth and remain spent at staging through seal. Staging anchors are (−15/2,−4,7173/2) and (−17/2,−17/4,7173/2) in comparison coordinates.

There are 69 actual local packet transfers using a 64-bit packet: namespace16 (0x4e41), operation-class8, body-index8, generation8, hop-ordinal16, completed1, inverse-address1 and six zero reserves. A valid copy requires the named actual incoming reciprocal port, predecessor, code and ordinal; wrong ports, absent predecessors, stale ordinals, noncode and halt cannot advance a valid chain.

Three disjoint retained 64-bit baton banks prevent mailbox reuse from erasing history. Across transfers, 4,416 distinct receiver words retain their packet/ordinal through seal. Local reciprocal ACK clears only reversible comparison scratch while its inputs remain; no remote reverse-cleanup oracle is assumed. A failed ACK/clear inhibits the dependent receipt.

Each packet transaction has a 1,024-scalar-slot microsequence: code checks, copies, equality/conjunction, ordinal/carry, inverse/ACK and guards. Each slot has its distinct local heat pack. Courier hops pay 128 custody posts, 512 scalar guard slots and 513 complete local heat packs. Literal paths, owners, port IDs, program/fault loading, schedule and retained histories are in the ACK data, not inferred from endpoint proximity.

## 9. UNA-65 atomic OFFER and independent sealing

### 9.1 Why an atomic counted event is specified

Reading 65 changing bits in sequence does not measure one instantaneous count across carry. Holding the independent clock until an apparatus is ready would defeat the cutoff. The design therefore adopts the fixed finite native operation **UNA-65**, an explicit additional branch premise, rather than either behavior.

The actual complete operation has **207 scalar operands / 828 dimers**, plus a distinct 38-dimer positive bank. The 140-word simultaneous payload label is not its complete arity. [acquisition-una65.json](data/acquisition-una65.json) lists every source, cutoff, latch, readiness, handshake and output operand.

Important address ranges are:

| Operand | Physical observer dimers / source |
|---|---|
| 65-bit independent counter | 0…259, four dimers per bit |
| 65-bit immutable cutoff | 260…519 |
| OPEN | 660…663 |
| ARRIVED | 664…667 |
| SEALED | 684…687 |
| New 65-bit latch | 3042…3301 |
| New latch-present | 3302…3305 |
| New preoffer-ready | 7914…7917 |
| New offer-event source | 7918…7921 |
| Head FINAL/OFFERED/inverse | Actual guide head dga(17920,0) ready-primary dimers 15520…15531 |

### 9.2 Complete bijection, guards and cost

For every code value n,x in 0…2^65−1, the operation exchanges disjoint matched blocks:

~~~
(n, latch=x, ARRIVED=0, present=0,
 head FINAL/OFFERED/inverse=1/0/0,
 offer-event=0, event-output=00)
↔
(n, latch=x XOR n, ARRIVED=1, present=1,
 head FINAL/OFFERED/inverse=0/1/1,
 offer-event=1, event-output=01).
~~~

Read-only readiness, OPEN and SEALED guards must have their specified valid values. All unmatched code sectors and every noncode/failed-readiness sector have the full identity extension on the capture payload; the separate independent clock still advances. The involution S satisfies S=S†=S^-1, so K=(I−S)/2 is positive, norm 1 and scalar/AM-invariant.

Forward dispatch additionally requires zero latch, no previous offer, exact final-ready completion, correct held-bay/port guards and current counter **strictly below cutoff**. The reciprocal is available only at the retained inverse address, not as a forward retry. The count source is nondemolition. Later counter increments never redispatch S or overwrite the captured latch.

The dedicated positive bank has 38 dimers, inlet 55k0, full ON/halted capacity 174k0; it is not borrowed from a nominal spare reserve. One π service error and 32 dispatch/guard posts are charged. These bookkeeping slots do not turn the indivisible admitted native capture into a continuous interruptible implementation. Uniform capture deviation is bounded by πζ_service on the full retained domain, while preserving the separate counter projectors/cadence.

### 9.3 Priority, failed attempts and one seal

At a simultaneous cutoff/OFFER, cutoff closes OPEN first; the tie is rejected. The clock is never held by preparation, acquisition or failure. Accepted capture records the one instantaneous native count in all lanes, not a file timestamp or an interpolated time.

After cutoff, the observer executes its prescribed acquisition suffix whether or not an offer arrived. It copies the original 38 snapshots and 512 descriptors once, then 512 source ternary diagnoses and 59 ordinary event snapshots. Events 59…62 write their actual native-transition outputs; reserved63 stays zero. There is exactly **one** final seal at the final cursor.

A failure midway through acquisition leaves event61 zero, but the independent clock still seals at its fixed final cursor. Event62 and SEALED can therefore be one in an unavailable record. Event61 is a success comparison, **not the seal-firing guard**. Invalid/noncode new A diagnostics cannot hold the C clock; their code-X takes its complement while the clock seal proceeds. A late completion or offer is retained separately and cannot replace the original attempted row.

### 9.4 Independent UNA-65 failure conditions

The premise fails if any of these occur in its intended domain: lane skew across a carry; clock hold caused by preparation/offer/failure; accepted cutoff tie; second accepted offer; subsequent latch overwrite; wrong cursor source port; noncode accepted as valid; unaccounted energy/AM transfer; or inability to realize the fixed bounded joint native completion. Producing a desirable number does not validate the operation.

## 10. Complete acquired alphabet and strict record schema

### 10.1 Output census

Each scalar occurrence is an observation of a named four-dimer output word. The complete ordered output mapping is [output-order.json.gz](data/output-order.json.gz); physical inputs remain a separate ordered list of 170 factors and four actual held heads.

| Field | Scalar occurrences | Full output dimers |
|---|---:|---:|
| Original 38 ternary snapshots | 76 | 304 |
| W/R three bits plus present | 4 | 16 |
| Four 128-bit descriptor destinations | 512 | 2,048 |
| Four×128 ternary source diagnoses | 1,024 | 4,096 |
| 64 ternary event diagnoses | 128 | 512 |
| Atomic offered count | 65 | 260 |
| Offer-count present | 1 | 4 |
| Literal SEALED | 1 | 4 |
| **Total** | **1,811** | **7,244 distinct dimers** |

Retained phase8 work is not an added phase sensor. No caller-supplied lifecycle Boolean, expected-ID table, archive filename or final cursor supplies a missing observation.

### 10.2 The original 38 words

The exact order is:

1. source_setting0, source_setting1, source_setting2;
2. receiver_setting0, receiver_setting1, READ;
3. PRESENCE, CUT, ROUNDTRIP;
4. endpoint.S_t, S_q, T_t, T_q, D_t, D_q, A_a, A_b, B_a, B_b;
5. endpoint.blank, quadrature, direction, n0,n1,n2,n3,n4, mode0,mode1,mode2;
6. endpoint.completed0,…,completed5, STOP0,STOP1.

Each word has two observed scalar bits: 00=u, 01=v, 10=retained source complement, 11=invalid. Each individual scalar itself has observed alphabet integer 0, integer 1, or literal "N" for noncode; N is not coerced to zero.

### 10.3 W/R and endpoint categories

W/R bits are read in their **printed order**, not reinterpreted with an integer endianness. With present=1, valid words map 001→positive, 010→negative, 011→none, 100→omit_READ; 101 is complement and 000 with present=0 denotes absence, neither an available angular result.

The ten endpoint flags, in order (S_t,S_q,T_t,T_q,D_t,D_q,A_a,A_b,B_a,B_b), define:

| Endpoint category | Ten-bit pattern |
|---|---|
| L | 1 0 0 0 0 0 0 1 1 0 |
| R | 0 0 1 0 0 0 1 0 0 1 |
| D | 0 0 0 0 1 0 0 1 1 0 |
| EMPTY | 0 0 0 0 0 0 0 1 1 0 |
| other | Any other **code-valid** ten-bit pattern |

There are 20 available joint categories: four angular outcomes × five endpoint categories, plus one unavailable/other category. “other” is a retained acquired pattern and forbidden; it is not silently thrown away. D is the local blank containing the arrived source; EMPTY means the endpoint received none, and its original source energy remains in the appropriate upstream control account.

### 10.4 Strict JSONL envelope

A delivered row has **exactly** the keys:

~~~
v, slot, kit, attempt, generation, case, provider,
input_refs, head_refs, observands
~~~

For intended slot i: v=1, kit=i, attempt=0, generation=0, case=i mod 821. Integers must be actual JSON integers, not booleans or floats. Provider has exactly version,ref, with version UNA65-1 and reference [i,0,0,"observer"]. Input/head arrays contain exactly 170 and four ordered references [i,0,0,physical_id].

Each scalar occurrence has exactly ref,value, ref [i,0,0,[dimer0,dimer1,dimer2,dimer3]] and value 0, 1 or N. Observands has exactly original_words, W_R, descriptor_destinations, descriptor_source_snapshots, event_snapshots, offer_count, offer_count_present, seal. Singleton fields are single occurrences, not one-element lists. All references are compared with the fixed native registry.

The envelope names the declared interface; it neither establishes a physical digitizer nor proves fresh independent kits. Physical occurrence provenance is a separate obligation, and the original archive is retained.

### 10.5 Classification order and all-attempt denominator

1. Begin every one of 65,680,000 intended slots as unavailable and forbidden.
2. Enforce exact envelope, scalar alphabets, physical reference/order/provider identity and one delivery per slot. Duplicate/reused/off-roster deliveries cannot win by last-record selection.
3. Decode all 38 ternary words and W/R. Any missing/malformed/noncode/incompatible output or reference remains unavailable+forbidden.
4. Compare all 512 token source/destination values and all 128 bits of each token, including all reserve bits and actual cargo counts.
5. Decode actual n from acquired endpoint words. Require tag events0…30 to match i<actual n, events31…62 to be one, event63 zero. Require present=1, valid atomic count below strict cutoff, SEALED=1, completed acquisition and retained custody.
6. If these independent acquisition checks pass, retain the acquired angular/endpoint category. Compare actual original/guide/endpoint controls, actual STOP and completed count with the prepared case; a wrong but valid value also marks the combined forbidden event.
7. Mark any category whose ideal support is zero, including code-valid endpoint “other,” forbidden.

Crucially, wrong-valid STOP or completed count **can retain an acquired category plus forbidden** when acquisition itself passes. Actual n versus observed tag-pattern inconsistency is acquisition-invalid. Coherently wrong n, count and tags relative to the scheduled case remains category+forbidden. Requested controls never fill missing observations.

Missing, partial, stopped, late, ties, warm/fault, malformed, noncode, duplicate, orphan and reuse records stay accounted. There is no successful-offer denominator or favorable retry. Nonempty structural_issues is a **separate admissibility failure**, even if the descriptive counts fit the numerical bands. Original archives and issue fields must remain attached to any result.


## 11. One complete error account

### 11.1 Fixed qualifications and final values

Use one good-map telescope over all source manufacture, angular controls, guide, material/cold/heat, endpoint, controllers, return, acquisition and native additions. The final coherent norm and nominal hot union are

ε ≤ 0.0961797078150192322826746794044…,

q_hot ≤ 1.34904794450741465030673538983×10^-7.

Service and bulk domain bounds remain 10^-15 and 10^-8. The final native addition contributes 1,646,728π service angle, 65,830,064 preoffer posts, 57,248 observer suffix posts and all its new stock; it does not consume an unused reserve as free work or earn a new probability allowance.

The final three claims are distinct:

| Quantity | Formula | Final upward bound |
|---|---|---:|
| Specified ideal-forbidden attempted event | ε²+d_in+q_hot+q_halt+d_dec | **0.009370671101** |
| General complete classical-record TV | ε+d_in+q_hot+q_halt+d_dec | **0.096299842720** |
| Predecoder same-factor whole-density half trace | d_in+ε+2√(q_hot+q_halt) | **0.116203193747704796439…** |

Only an ideal-annihilating event permits squaring coherent amplitude error. The joint density input bound enters **linearly**, not inside an amplitude square. The whole-custody statement keeps the actual physical purifiers on the same factors and uses a square-root bad-component allowance; it is not a ket bound or a postdecoder quantum fidelity claim.

The linear nominal hot/halt union at the reduced-record boundary requires the passive tag purifiers to remain untouched by later operations. Retaining coherent purifiers and then assuming their amplitudes are merely probabilities would give a different claim. No new stage or success-selected subensemble resets this accounting.

### 11.2 Interpretation and limiting cases

The 0.9370671101% number applies to the defined union of ideal-zero outcomes and attempted-operation failures under all qualifications. It does **not** mean every category, arbitrary input channel or retained physical state is accurate to better than 1%.

Blank/CUT/ROUNDTRIP have exact ideal dark behavior under their named operations; physical deviations still inherit the full attempted-event account. Removing the first endpoint mixer yields 1/2 independently of tags, while the phase panel retains the deterministic zero-tag and transfer cases. The leaf perturbation is nonzero; setting it to zero would give a different guide without the stated transverse construction. Independent within-spine detuning lies outside exact matching; freely fitting a detuning function to observed fringes is not allowed.

## 12. Frozen prospective comparison

### 12.1 Exactly 821 prepared cases

[prepared-cases.json.gz](data/prepared-cases.json.gz) gives every exact case and expected control word:

- **770 PHASE cases:** the registered endpoint phase/control roster, cycling the 21 bright original angular controls by j mod 21.
- **48 GUIDE_CONTROL cases:** every one of the 24 original angular controls under CUT and ROUNDTRIP, each with the fixed quiet endpoint word.
- **3 ORIGINAL_BLANK cases:** original pre-launch blank against the three analysis bases, guide PASS, same quiet endpoint word.

This is not a full Cartesian angular×phase experiment. The exact 21-bit prepared case field names the initial case; it does not manufacture actual A controls from an ideal metadata label.

For the registered angular q vector and endpoint probabilities p, the ideal joint cell is q_y p_z. The admitted timing/custody coordinates have fixed support conditions but no invented distribution within a timing band. Source data and phase selection are fixed independently of outcomes.

### 12.2 Counts, concentration and exact test

Every case has 80,000 intended attempts. There are 20 available joint-cell indicators plus the combined forbidden event per case, giving 821×21=17,241 nonredundant indicators. Unavailable is already included in the forbidden event. With slack η=1/100, 2nη²=16.

The finite exact rational family-tail certificate is

2·17,241 / Σ_(j=0)^100 16^j/j! < 0.003881.

The evaluated value is approximately 0.0038804378946694926, below the declared 0.01 family level. The finite positive series is a lower bound for e^16; its reciprocal therefore gives a rigorous upper tail enclosure. The exact rational numerator and denominator are in complete-accounts data.

A case is numerically compatible only when its forbidden fraction is at most

0.009370671101+0.01 = **0.019370671101**,

and every category fraction is within its exact ideal interval expanded by

0.096299842720+0.01 = **0.106299842720**.

Then enforce the separate structural admissibility gate. A comparison result remains an uncertified descriptive reduction until physical occurrence, actual digitization, full decoder qualification and repetition premises have their receipts.

Concentration requires **independently qualified complete fresh kits**, or a separately justified per-history conditional record law. One-shot marginals, unique IDs or a fresh-looking filename do not establish that premise. Shared latent faults can invalidate concentration while leaving individual marginal bounds unchanged.

This is a compatibility test, not observed sub-1% performance certification: a 1.5% forbidden frequency can pass. Categorywise bands are not a full empirical-TV certificate: ideal quarters shifted to .35,.35,.15,.15 have TV 0.2 although each category shifts only 0.1. The target is nevertheless discriminating: zero-tag deterministic L in case0 differs by 0.5 in its joint cell from a phase-independent half/half alternative, well beyond the category band.

### 12.3 Executable identity and archive custody

The sealed original target has TARGET SHA-256 fc812b86f4e90d38b6a15c4c5284010b8d1d6f4f9e0c024b0178504133e60d3b and externally pinned entry SHA-256 b722c5ddd50994982cc523353b6923a55e1ad3f4991d00f6819b3b511e9e7c75.

The original entry verifies source bytes before compiling the real adapter, checks the actual source/entry again afterward, and makes a final archive-identity read immediately before returning. Archive identity includes resolved path/device/inode/size/hash; source identity is bytes. This is bounded integrity under trusted runtime/filesystem and no concurrent mutation, not a hostile-race, atomic-snapshot or self-authenticating-code theorem. Hashes do not turn synthetic records into physical provenance. The curated exhibit data and any explanatory classifier are derivative references, not replacements for that exact entry.

## 13. Full resources, capacities and retained output inventory

### 13.1 Whole-system bill of quantities

These are sufficient **upper accounts**, not necessary lower bounds, a fabrication quote, physical dissipation or measured power.

| Item | Per fresh attempt |
|---|---:|
| Complete physical dimers | **31,007,572,988,523** |
| Complete independent observer | **43,146** dimers |
| Exact independent pure C | 674 factors |
| Pre-OFFER native posts P* | 1,504,551,385,726,580,037 |
| All native posts, including observer | **25,577,373,557,351,935,515** |
| Final cursor | 24,072,822,171,625,355,478 |
| Actual input word factors | 170 |
| Held measurement heads | 4, full 107,666-dimer bodies |
| Acquired output dimers | 7,244 distinct |

The final additive acquisition bill is:

| New family | Dimers |
|---|---:|
| Observer extension | 40,104 |
| 44 ordinary contacts | 169,136 |
| Eight complete raw support packs | 48 |
| Two complete receipt couriers | 215,332 |
| 823,477 fresh ten-dimer heat packs | 8,234,770 |
| **Total addition** | **8,659,390** |

The observer's new physical ranges are ready primaries 3042…9185; scalar remainders 9186…42977; cadence 42978…42999; modes 43000…43107; UNA positive bank 43108…43145. All are A, not additional exact C.

### 13.2 Complete positive inlet and capacity

The per-kit positive inlet bound in k0 units is exactly

~~~
1240810646036703567801522499374518385920239088876587792841370232215 / 4.
~~~

The full capacity bound, including active and halted interaction sectors, in k0 units is

~~~
19820215732868705230648580000241451461210988412447063837033573535839 / 16.
~~~

These include positive supplies, supports, programs, purifiers, work, heat, native controls and every retained failed sector. A completed good-path net energy balance cannot cancel installed positive capacity. Negative-time/negative-energy services, signed reserve credits and discarded failed stores are not permitted.

For the guide-only edge/vertex assembly, conservative transport+record primitive stock is 2,734,080 dimers before full material/controllers/preparation, with 205,607,945,504 posts. That is a **subassembly** number, not the whole-device bill. Its smaller alternative capsule census is not imported into the final qualified construction. The same distinction applies to any source, instrument or endpoint subtotal.

### 13.3 Complete no-reuse campaign

| Quantity | 65,680,000 fresh attempts |
|---|---:|
| Complete physical dimers | **2,036,577,393,886,190,640,000** |
| Native posts | **1,679,921,895,246,875,124,625,200,000** |
| Independently pure C factors | 44,268,320,000 |
| Nonclock A dimers | 2,036,577,393,841,922,320,000 |

Exact campaign inlet and capacity are the per-kit rationals multiplied by 65,680,000 and are provided in [complete-accounts.json](data/complete-accounts.json). No time in seconds follows from these counts. Fresh-kit repetition is finite provisioned operation, not a promised indefinitely reusable physical fiber.

### 13.4 What remains after one attempt

The final physical output includes the original source in its returned berth with its actual state/reference, its source-side partner/copy/M stock, spent guide/endpoint work, all three angular components and recoil, measured pointers/phase workspace, four held measurement heads and cargo, distinct returned empty controllers, spent receipt couriers at staging, original/new receipt and baton banks, active or halted switches/tickets, all used/unused heat stock, formation/contact/outlet records, fault fields, purifiers and unreturned/failed bodies. A missing body is a failed custody outcome, not a vanished resource.

Good work-return or positive-switch return is not a full reset. Every uncomputation reverses only the named reversible scratch with retained inputs. It never deletes a measurement, history, purifier or consumed reference. The next attempt consumes a separately qualified kit unless a new complete reset/reuse law is derived and reviewed.

## 14. Material and SI engineering determination register

The finite laws and dimensionless choices above are assigned. Their actual material realization and useful engineering performance are not.

| Determination still required | Fixed role / independent qualification | What cannot substitute |
|---|---|---|
| Physical circulation/material family | Realize the full dimer/scalar/contact catalogue with its complete complements and retained conservation accounts | A familiar particle, glass or vendor label |
| Work scale k0 | Independent same-account finite work/translation-tariff calibration | Choosing k0 after the fringe to manufacture optical frequency |
| Length scale a0 | Independent reciprocal geometry/rest-frame comparisons | Reading drawing pixels as micrometers |
| Reference time / cadence | Independent actual native completion-count comparison across the whole horizon | File timestamps, final cursor, or offer-induced clock hold |
| LCR matching and isolation | Exact common-coordinate ownership relations, stationary interval and scalar shape/load behavior | Unmeasured small detuning fitted to output |
| Formed-body geometry | Complete finite coordinates, partner/port/rest relations and correlated body domain | Arbitrary independent microscopic jitter or continuous swept-volume assumptions |
| Cold joint input | One full A/reference half-trace bound before preparation with retained purifiers | Marginal coldness, a pure-looking output or success-selected runs |
| Thermal family | Fixed gap/pair kernel, local opportunity count, OFF windows and retained passive tags | An unspecified ambient bath or unpriced heat sink |
| Native completion and UNA-65 | Fixed finite local cadence, reciprocal custody, atomic count, priority and full error domain | A general universal-unitary assumption or sampled sequential bits |
| Actual digitization | Fixed map from actual occurrences to the exact alphabet, complete classical decoder TV, missing/noncode/partial handling | Names, hashes, expected output arrays or a provider-qualified Boolean |
| Repetition | Whole-kit independence or a justified per-history record law | Unique kit identifiers or one-shot marginal probability bounds |
| Practical resources and SI performance | New version-correct reduction of stock/posts and measured bandwidth/latency/energy/length | Multiplying a dimensionless phase by an assumed carrier frequency |

Thus no fiber diameter, jacket, refractive-index profile, attenuation in dB/km, operating wavelength, bandwidth, bit rate, receiver noise floor, clock oscillator model, detector electronics, temperature in kelvin, cryostat, laser, connector, laboratory supply or vendor part is assigned by this paper. Those are not blank fields to fill arbitrarily: each must be chosen against the fixed native observable and material law and then carry its own engineering and measurement qualification. The theoretical positive inlet/capacity is not a wall-plug energy or heat-load estimate.

## 15. Failure and acceptance requirements

The conditional construction is accepted as this version only when the following conjunction is preserved:

1. **Same source:** exact source/copy/M targets, complete prepared seed, local command ordering, same guide/receiver factor, full SWAP, actual ra24 endpoint and actual return. No replacement signal.
2. **Same physical domain:** full nine-state dimers, three angular components, exact support membership/rate law, complete positive local operators and noncode/wrong-loading complements.
3. **Same material and geometry:** formed correlated body templates, all actual contacts/ports/feet/rest vectors, scalar isolation, actual finite vacancy/custody routes and local service radius. No zero-cost reshaping, compressed controller or unproved transport between apparently nearby points.
4. **Same complete accounts:** source, terminal, both work pairs, retained transverse/support energy, positive ON/halted stores, heat, reactions, programs, return/reset/failure and every purifier.
5. **Same control/clock:** named finite native catalogue and fields, one P* horizon/fault denominator, bounded progress, local packets and actual acknowledgments, exactly one atomic counted OFFER and unconditional final seal.
6. **Same record:** 170 original factors, 38 words, four held heads, 1,811 scalar outputs/7,244 dimers, exact token/reserve/event relations, correct D versus EMPTY and all intended attempts retained.
7. **Same error definitions:** one initial joint half-trace, exactly 674 independent pure C factors, full-domain service/bulk bounds, one halt and one complete decoder bound; event, record-TV and predecoder density not conflated.
8. **Same comparison:** 821-case roster, 80,000 attempts per case, fixed simultaneous indicators/slack, requisite repetition premise, separate structural admissibility, original archives and exact externally pinned source-bound entry.

Every detailed native-author/crosswalk, physical-review, accounts-review, raw-author/review, inherited adopted and final-goal condition remains conjunctive; this short register cannot waive a source condition. Changing an operation, physical domain, stock, geometry, timing, records or claimed performance creates a new version needing its corresponding proof and review.

### Engineering consequence

The construction shows how one finite, explicitly supplied circulation organization could carry relational angular information over a guided route, expose a paid phase-dependent endpoint result, and retain a causally complete failure-inclusive record without replacing its source or losing its energy and angular-momentum accounts. It gives a concrete mathematical target for material and resource engineering. The next useful engineering advance is not a stronger label for this result, but a realizable qualified type and a version-correct reduction of its very large supply and execution requirements.

## 16. Exact engineering companions

The downloadable derivative engineering archive includes:

- all 2,049 guide ROM integers and guide native field/mode/CUT bindings;
- all 821 prepared cases with original/guide/endpoint control words and ideal category laws;
- all 170 input references, four held-head references, 1,811 ordered scalar output references and complete raw schema;
- all 207 named UNA-65 scalar operands and complete matched/unmatched/complement rules;
- all four exact 128-bit token bodies, source/destination/snapshot references and initial/final bit patterns;
- all 64 event sources, outputs, writers and completion dependencies;
- actual receipt-courier paths, 69 transfer sequences, 128 courier-hop schedules and physical ownership;
- source, endpoint-route, root and acquisition contact/rest data, native field lifetimes and heat-atlas ownership;
- exact final per-kit/campaign resource rationals, errors and statistical constants;
- dependency-linked current geometry, complete retained annex BOM/endpoint generators, preparation/compiler and final-classifier definitions, explicit module bindings, exact retained support tables and a section/hash provenance map.

All data companions are explicitly **derivative reference material**, not a new sealed runtime or physical data. The original source collection is unchanged. Arithmetic and mapping checks for this appendix are recorded in [DATA-VERIFICATION.json](DATA-VERIFICATION.json); they verify the exposition and byte-derived data, not physical material, hardware or observations.


### 16.1 Dependency and exact-number reading contract

[normative-functions.json.gz](data/normative-functions.json.gz) contains the retained current source definitions, including the optical annex `stats`, `roles`, `caps`, `ordinary_contacts` and `atlas_point` geometry/BOM family, their local-factory geometric utilities, the guide spine feet, endpoint ledgers, acquisition model/relay/atlas and the final raw classifier. [NORMATIVE-DEPENDENCIES.json](data/NORMATIVE-DEPENDENCIES.json) resolves local imports and dynamic aliases; [normative-support-data.json.gz](data/normative-support-data.json.gz) supplies exact retained values for the source selectors those definitions use. A source hash by itself is not an engineering definition. The module and literal-value bindings provide the definitions here; the hashes establish their provenance. Only named integrity/publication/diagnostic or superseded semantic symbols are excluded, with their reasons listed. In particular, `spine_foot_ref` and `spine_foot` are physical definitions and are retained.

These are readable engineering references, not a newly sealed executable. They describe component functions in their declared namespace; use the final complete-account tables, actual 821-case roster and final adapter whenever a component's inherited wrapper or earlier subtotal differs. The original 72/770-case semantic classifiers and historical target wrappers are not the final acquisition protocol. No private input manifest or historical repair bundle is needed to interpret the delivered current definitions.

**Exact integer handling is mandatory.** Native cursor, token, resource and campaign values exceed 2^53−1. Read numeric JSON integers with an arbitrary-precision parser (for example Python's ordinary `json.load`, whose integer values are Python `int`), and rational/decimal strings with exact arithmetic. Ordinary JavaScript `JSON.parse` converts numeric tokens to `Number` and can silently round them before a reviver runs; use a lossless parser instead. A web-distribution derivative may encode large integers as decimal strings; convert those strings directly to arbitrary-precision integers, never through `Number`. Preserve original source bytes and source hashes separately from any derivative encoding. This table-reader rule does **not** relax or change the final raw-record grammar: fields required to be actual JSON integers remain actual integers in raw acquisition records.
