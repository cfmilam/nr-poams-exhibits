#!/usr/bin/env python3
"""
The Valley — a conditional audit of the smooth single-beta energetic floor at
fixed A using inherited model coefficients. Companion script to valley-of-stability.html
(NR/POAMS exhibit series). Run 2026-09-06; integer-score repair 2026-09-08.

THE QUESTION. At fixed census A, the quanta split between the pattern channel
(count Z; extant label: protons) and the core twist ledger (count A-Z; extant:
neutrons). The model reads a one-step beta move as a re-filing between the
channels at fixed A. This instrument evaluates a static cost surface. It does
not calculate a weak decay rate, lifetime, or stability against other modes.
Its coefficient packages predate this program, but the 2026-09-08 Mass Ledger
audit quarantines the cluster calibration behind the imbalance magnitude.
Accordingly the numerical result is conditional, not an independent
coefficient derivation:

  CENSUS-IMBALANCE COST (Mass Ledger P3 + section 6 candidate):
    form (A-2Z)^2/A structurally motivated by the two filing ladders and
    sense-relabel evenness. The displayed coefficient package is
    a_sym = (1/3)E_F + C/(1+2*gamma) = 24.0-24.9 MeV
      - (1/3)E_F = 11.13 MeV: filling-ladder tax at capacity density
        rho0 = 3/(4 pi r0^3), a model-selected capacity density
      - C/(1+2g): contact column frozen at capacity; C = a_v + tau_b = 35.85
        (quarantined cluster calibration); gamma = 0.85(5) (named import)
  FAR-FIELD CIRCULATION COST (Mass Ledger P5 candidate):
    a_c * Z(Z-1)/A^(1/3), a_c = (3/5) * k / r0 = 0.72 MeV within the
    uniform-sphere approximation and rounded inputs. Z(Z-1) is pair counting;
    k = 1.44 MeV fm, the strain quantum (= alpha*hbar*c -- alpha enters ONLY
    as the corpus's declared contingent input); r0 = 1.2 fm (named import).
    Booked exchange filing layer (T-P3', untuned): -0.53 * Z^(4/3)/A^(1/3).
  PAIRING (Mass Ledger P4): parity term, ee < eo < oo forced; magnitude class
    Delta ~ 12/sqrt(A) used ONLY as tagged comparative shadow for byproduct
    thresholds (the ledger's own plateau count is its named open item 1).
  CANDIDATE DEVELOPED HERE (F1b): the confined-ladder re-filing seat for the kinetic
    part of a_sym -- two classes fill separate copies of the Inner Census
    ladder N = 2 n_r + l, capacities (N+1)(N+2); continuum E = (h Dnu/4)(3n)^(4/3)
    per class gives a_kin = h*Dnu*A^(1/3) / 144^(1/3), and with the inherited
    spacing h*Dnu ~ 41.7 A^(-1/3) MeV (eyewall_pricing.py byproduct)
    the A^(1/3) cancels: a_kin(ladder) = 7.96 MeV, A-independent.
    Cross-seat tension NAMED: 7.96 (ladder) vs 11.14 (capacity density), the
    oscillator's skin-spread top rungs vs the pinned interior density.

ORDER NOTE (honesty): the original registration and the 2026-09-08 repair
registration print before computation. The repair was prompted by a known
referee result: the old gate compared fractional minima with integer nuclides.
It is therefore an audit, not a new blind prediction. No coefficient is tuned
on this page; several coefficients are measured or model-selected imports and
remain tagged as such.

UNITS (corpus convention, ratified 2026-08-15): displays print turn rate nu
on the Hz ladder (1 record MeV = 241.8 EHz = 0.2418 ZHz of whole closures).
Internal arithmetic runs in the record's MeV because every import is quoted
there; every PRINTED coefficient leads native with the record value tagged.
"""
import math

EHZ_PER_MEV = 241.799
def zhz(x): return x * EHZ_PER_MEV / 1000.0
def ehz(x): return x * EHZ_PER_MEV

# ---------------- named ledger constants (all cited above) ----------------
HBAR2_OVER_M = 41.5        # MeV fm^2
R0 = 1.2                   # fm
AV = 15.75                 # MeV
CCAP = 35.85               # MeV (quarantined inherited cluster calibration)
K_STRAIN = 1.44            # MeV fm (strain quantum = alpha*hbar*c; contingent input)
GAMMA = (0.80, 0.85, 0.90) # named import band
A_C = 0.6 * K_STRAIN / R0  # = 0.72 MeV within the uniform-sphere model
EXCH = 0.53                # MeV, booked exchange filing layer (T-P3', untuned)
CLASS_STEP = 0.782         # MeV: n heavier than p+e (measured record import,
                           # TAGGED: the class rest split is NOT yet corpus-derived)
DELTA_PAIR = 12.0          # MeV/sqrt(A): pairing shadow, tagged comparative only
RHO0 = 3.0 / (4.0 * math.pi * R0**3)

# ---------------- registration (printed before any valley number) ----------
print("=" * 78)
print("THE VALLEY -- INTEGER-SCORE REPAIR (registration first, comparison last)")
print("=" * 78)
print("""REPAIRED GATES (2026-09-08; target-known audit, not preregistration):
 G1  minimize over INTEGER Z and require distance <= 2 charge units from the
     nearest selected evaluated one-step-beta reference at
     A = {16,40,56,120,208,238}. Reference sets are
     {8}, {18,20}, {26}, {50,52}, {82}, {92}. A=40 and A=120 each carry two
     even-even local one-step-beta minima; U-238 is beta-stable but alpha
     radioactive. +-2 is MODEL TOLERANCE, not uncertainty in integer Z.
 G2  the CONTINUOUS ENVELOPE bends: Z_env/A falls from ~0.48 (A=16) to
     0.37-0.40 (A=238) in EVERY package (record: 92/238 = 0.3866). The bend
     scale is the conditional output of the displayed ratio a_c/4a_sym; no
     coefficient is adjusted within this script.
 G3  MODEL BYPRODUCT (report only): the tagged pairing shadow can create two
     even-even local one-step-beta minima separated by an odd-odd bridge.
     This is energetic permission only: no rates, lifetimes, branching ratios,
     or stability against alpha/fission/double-beta decay are computed.
FORK TABLE (both a_sym seats are inherited model candidates; the strain
dressing terms are booked/tagged; the four packages are a consistency cross):
 P-I   a_sym = frozen-capacity seat (5.79-6.03 ZHz; record 24.0-24.9 MeV)
       + bare P5 strain [historical conditional pairing of terms]
 P-II  a_sym = ladder re-filing seat (5.02-5.26 ZHz; record 20.8-21.7 MeV)
       + dressed strain
       (P5 + exchange layer + class rest step)                   [full-dress]
 P-III a_sym = frozen-capacity + dressed strain   (mixed corner)
 P-IV  a_sym = ladder re-filing + bare strain     (mixed corner)
KNOWN REPAIR EXPECTATION: integer scoring makes all four packages pass the
 inclusive +-2 gate. Therefore G1 cannot discriminate the package pairing.
""")

# ---------------- coefficient assembly (cited, no freedom) ----------------
def e_fermi_third():
    kf = (3.0 * math.pi**2 * RHO0 / 2.0) ** (1.0 / 3.0)
    ef = HBAR2_OVER_M * kf * kf / 2.0
    return ef / 3.0

def contact_band():
    return tuple(CCAP / (1.0 + 2.0 * g) for g in GAMMA)  # decreasing in gamma

# confined ladder (Inner Census): exact spacing byproduct + re-filing seat
def rung_cap(N): return (N + 1) * (N + 2)

def ladder_fill_weight(nq):
    tot, N, w = 0, 0, 0.0
    while tot < nq:
        take = min(rung_cap(N), nq - tot)
        w += take * (N + 1.5)
        tot += take
        N += 1
    return w

def hw_of_A(A):
    nq = A // 2
    msr = 0.6 * (R0 * A ** (1 / 3.0)) ** 2
    b2 = msr * nq / ladder_fill_weight(nq)
    return HBAR2_OVER_M / b2

def ladder_energy(nq, hw):
    """Exact discrete ladder energy of nq quanta of one class (units MeV)."""
    tot, N, e = 0, 0, 0.0
    while tot < nq:
        take = min(rung_cap(N), nq - tot)
        e += take * (N + 1.5) * hw
        tot += take
        N += 1
    return e

def ladder_akin_exact(A):
    """Discrete two-ladder re-filing coefficient at census A. NOTE (booked as
    structure, not a bug): a pure rung is degenerate, so re-filing WITHIN a
    partially filled rung is locally free -- the quadratic cost lives in the
    seam staircase. The coefficient below is the seam-averaged curvature:
    least-squares fit of E(x)+E(-x) to c0 + c2*(2x)^2/A over a window wide
    enough to span several rung seams (x = (N-Z)/2 in [-48, 48])."""
    hw = hw_of_A(A)
    n0 = A // 2
    xs = [x for x in range(-48, 49, 2) if 0 < n0 + x < A]
    es = [ladder_energy(n0 + x, hw) + ladder_energy(A - (n0 + x), hw) for x in xs]
    qs = [(2 * x) ** 2 / A for x in xs]
    qb = sum(qs) / len(qs)
    eb = sum(es) / len(es)
    num = sum((q - qb) * (e - eb) for q, e in zip(qs, es))
    den = sum((q - qb) ** 2 for q in qs)
    return num / den, hw

AKIN_CAP = e_fermi_third()
C_LO, C_MID, C_HI = sorted(contact_band())
HW41 = hw_of_A(238) * 238 ** (1 / 3.0)          # the byproduct constant ~41.2
AKIN_LAD = HW41 / (144.0 ** (1.0 / 3.0))        # continuum: A-independent

print("conditional coefficient assembly (inherited inputs; none adjusted here;")
print("native turn rate first, record MeV tagged):")
print(f"  rho0 = {RHO0:.5f} fm^-3;  (1/3)E_F (capacity-density seat) = "
      f"{zhz(AKIN_CAP):.3f} ZHz (record {AKIN_CAP:.2f} MeV)")
print(f"  contact frozen at capacity C/(1+2g), g={GAMMA}: "
      f"{zhz(C_LO):.2f} / {zhz(C_MID):.2f} / {zhz(C_HI):.2f} ZHz "
      f"(record {C_LO:.2f} / {C_MID:.2f} / {C_HI:.2f} MeV)")
print(f"  ladder spacing byproduct h*Dnu*A^(1/3) = {zhz(HW41):.2f} ZHz "
      f"(record {HW41:.1f} MeV; booked ~41)")
print(f"  ladder re-filing seat (continuum): a_kin = {zhz(HW41):.2f}/144^(1/3) "
      f"= {zhz(AKIN_LAD):.3f} ZHz (record {AKIN_LAD:.2f} MeV), A-independent "
      f"(the A^(1/3) cancels)")
for A in (56, 120, 238):
    ak, hw = ladder_akin_exact(A)
    print(f"    exact discrete check A={A}: a_kin = {zhz(ak):.2f} ZHz "
          f"(record {ak:.2f} MeV; h*Dnu = {zhz(hw)*1000:.0f} EHz)")
ASYM_CAP = (AKIN_CAP + C_LO, AKIN_CAP + C_MID, AKIN_CAP + C_HI)
ASYM_LAD = (AKIN_LAD + C_LO, AKIN_LAD + C_MID, AKIN_LAD + C_HI)
print(f"  a_sym  frozen-capacity seat: {zhz(ASYM_CAP[0]):.2f} .. {zhz(ASYM_CAP[2]):.2f} ZHz "
      f"(central {zhz(ASYM_CAP[1]):.2f}; record {ASYM_CAP[0]:.2f}..{ASYM_CAP[2]:.2f} MeV)"
      f"   [Mass Ledger section 6; quarantined cluster calibration]")
print(f"  a_sym  ladder seat:          {zhz(ASYM_LAD[0]):.2f} .. {zhz(ASYM_LAD[2]):.2f} ZHz "
      f"(central {zhz(ASYM_LAD[1]):.2f}; record {ASYM_LAD[0]:.2f}..{ASYM_LAD[2]:.2f} MeV)")
print(f"  a_c = (3/5)k/r0 = {ehz(A_C):.1f} EHz (record {A_C:.3f} MeV; "
      "exact only within the displayed uniform-sphere model and rounded inputs);")
print(f"  exchange layer {ehz(EXCH):.1f} EHz (record {EXCH} MeV) * Z^(4/3)/A^(1/3);")
print(f"  class rest step {ehz(CLASS_STEP):.1f} EHz (record {CLASS_STEP} MeV) "
      f"(tagged measured import)")
print()

# ---------------- the valley: cost function and minimizer -----------------
def cost(Z, A, asym, dressed):
    E = asym * (A - 2 * Z) ** 2 / A + A_C * Z * (Z - 1) / A ** (1 / 3.0)
    if dressed:
        E -= EXCH * Z ** (4.0 / 3.0) / A ** (1 / 3.0)
        E -= CLASS_STEP * Z          # re-filing N->Z returns the class step
    return E

def zstar_cont(A, asym, dressed):
    lo, hi = 1.0, A - 1.0
    for _ in range(200):
        m = (lo + hi) / 2
        d = cost(m + 1e-4, A, asym, dressed) - cost(m - 1e-4, A, asym, dressed)
        if d < 0: lo = m
        else: hi = m
    return (lo + hi) / 2

def zstar_int(A, asym, dressed):
    return min(range(1, A), key=lambda Z: cost(Z, A, asym, dressed))

PACKAGES = [
    ("P-I   frozen-capacity + bare ", ASYM_CAP, False),
    ("P-II  ladder + dressed       ", ASYM_LAD, True),
    ("P-III frozen-capacity+dressed", ASYM_CAP, True),
    ("P-IV  ladder + bare          ", ASYM_LAD, False),
]
ANCHORS = [16, 40, 56, 120, 208, 238]

print("MODEL COMPUTATION -- Z*(A) per package (continuum; [band from gamma];")
print("integer minimum in parens). The comparison target is known; record prints later.")
rows = {}
for name, band, dressed in PACKAGES:
    out = []
    for A in ANCHORS:
        zc = zstar_cont(A, band[1], dressed)
        zl = zstar_cont(A, band[0], dressed)
        zh = zstar_cont(A, band[2], dressed)
        zi = zstar_int(A, band[1], dressed)
        out.append((A, zc, min(zl, zh), max(zl, zh), zi))
    rows[name] = out
    cells = "  ".join(f"A={A}:{zc:5.2f}[{a:.1f},{b:.1f}]({zi})" for A, zc, a, b, zi in out)
    print(f" {name}: {cells}")
print()
print("bend check: Z*/A across A (central values):")
for name, band, dressed in PACKAGES:
    line = "  ".join(f"{A}:{zstar_cont(A, band[1], dressed)/A:.3f}" for A in ANCHORS)
    print(f" {name}: {line}")
print()

# byproduct G3: parity staircase thresholds (pairing shadow, tagged)
def curv(A, asym): return 8 * asym / A + 2 * A_C / A ** (1 / 3.0)
print("G3 model byproduct (pairing shadow Delta = 2.90/sqrt(A) ZHz -- record")
print("12/sqrt(A) MeV -- tagged comparative):")
oo_thresh = next(A for A in range(10, 200)
                 if curv(A, ASYM_CAP[1]) / 2.0 < 2 * DELTA_PAIR / math.sqrt(A))
print(f"  model odd-odd local-minimum window closes at A ~ {oo_thresh} "
      f"(condition: half-curvature step < 2*Delta)")
def one_step_local_minima(A, asym=ASYM_CAP[1], dressed=False):
    """Local integer minima on the smooth parabola + tagged parity shadow.
    This tests energetic permission for one-step beta moves only; it predicts
    neither decay rates nor stability against other channels."""
    D = DELTA_PAIR / math.sqrt(A)
    zc = zstar_cont(A, asym, dressed)
    def par(Z):
        N = A - Z
        if A % 2 == 1: return 0.0
        return -D if Z % 2 == 0 else +D
    Zs = range(max(1, int(zc) - 5), int(zc) + 6)
    E = {Z: cost(Z, A, asym, dressed) + par(Z) for Z in Zs}
    return [Z for Z in list(Zs)[1:-1] if E[Z] < E[Z - 1] and E[Z] < E[Z + 1]]

first_doublet = next(A for A in range(10, 300, 2) if len(one_step_local_minima(A)) >= 2)
bands = [(20, 60), (60, 120), (120, 200), (200, 260)]
frac = {b: sum(1 for A in range(b[0], b[1], 2) if len(one_step_local_minima(A)) >= 2)
            / max(1, len(range(b[0], b[1], 2))) for b in bands}
print(f"  first even-A pair of local one-step minima: A ~ {first_doublet}; model")
print("  model fraction by band: " + ", ".join(f"A {a}-{b}: {frac[(a,b)]:.0%}" for a, b in bands))
print("  odd A: the parity term is class-blind, so integer sampling normally")
print("  gives one energetic floor (two only at an exact half-integer tie).")
print("  No observational-stability or lifetime claim is made from this output.")
print()

# ---------------- the gate, applied last ----------------------------------
print("=" * 78)
print("COMPARISON (the gate, applied last)")
print("=" * 78)
Z_REF = {16: (8,), 40: (18, 20), 56: (26,), 120: (50, 52),
         208: (82,), 238: (92,)}
print("selected evaluated one-step-beta reference charges (tagged): "
      + ", ".join(f"A={A}: Z={refs}" for A, refs in Z_REF.items()))
print("  U-238 is included as beta-stable but is alpha radioactive; these six")
print("  chains are a sparse orientation set, not a full stability census.")
verdicts = {}
for name, band, dressed in PACKAGES:
    devs = []
    for A, zc, zl, zh, zi in rows[name]:
        refs = Z_REF[A]
        nearest = min(refs, key=lambda z: abs(zi - z))
        devs.append((A, zi - nearest, zi, nearest))
    mx = max(abs(d) for _, d, _, _ in devs)
    ok = mx <= 2.0
    verdicts[name] = (ok, mx)
    cells = "  ".join(f"A={A}:{zi}->{ref} ({d:+d})" for A, d, zi, ref in devs)
    print(f" {name}: {cells}  max|dZ|={mx:d}  G1 {'PASS' if ok else 'FAIL'}")
print()
bend_ok = all(0.37 <= zstar_cont(238, band[1], dressed) / 238 <= 0.40
              and zstar_cont(16, band[1], dressed) / 16 > 0.45
              for _, band, dressed in PACKAGES)
print(f"G2 bend (all packages): continuous Z_env/A(238) in [0.37,0.40] "
      f"and Z_env/A(16) > 0.45: {'PASS' if bend_ok else 'FAIL'}")
print(f"G2 record value: 92/238 = {92/238:.4f}")
print()
print("VERDICT:")
print(" - The bend is package-independent within this audit: it follows from")
print("   a_c/4a_sym conditional on the displayed booked/imported coefficients.")
print(" - Integer G1 does NOT discriminate the four packages; all pass:")
for name in rows:
    ok, mx = verdicts[name]
    print(f"     {name}: {'PASS' if ok else 'FAIL'} (max {mx:.2f})")
print(f" - The earlier continuous-score package discriminator is retired. The")
print(f"   kinetic a_sym seat (capacity-density "
      f"{zhz(AKIN_CAP):.2f} ZHz")
print(f"   vs ladder re-filing {zhz(AKIN_LAD):.2f} ZHz; record {AKIN_CAP:.2f} vs "
      f"{AKIN_LAD:.2f} MeV) is correlated with the strain dressing;")
print("   the sparse smooth-valley gate cannot select their package. The")
print("   record-implied composite (~5.66 ZHz bare / ~5.32 ZHz dressed;")
print("   record ~23.4 / ~22.0 MeV) is inferred from the same target, not held-out evidence.")
print(f" - NAMED IMPORT: the class rest step {ehz(CLASS_STEP):.1f} EHz "
      f"(record {CLASS_STEP} MeV; measured; un-derived).")
