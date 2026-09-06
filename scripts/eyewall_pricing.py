#!/usr/bin/env python3
"""
Eyewall pricing — the Inner Census's declared open joint, executed 2026-09-06.

Companion script to inner-census.html (NR/POAMS exhibit series). FROZEN as the
judge of the intruder threshold: the gate arithmetic and the l* = 3 band below
must never be altered; new candidate mechanisms may only be ADDED as separate
pricing functions and run through the same blind gate.

VERDICT (booked): FALSIFICATION — under the corpus-licensed coupling class the
native pricing yields no threshold near l* = 3 (miss ~45x at l = 3, bare sign
anti-record). Per the exhibit's own falsifier, the intruder counting rule
stands demoted to a description of the record. See inner-census.html §6.

Prices the co/counter rate split kappa(l) against the confined ladder's rung
spacing, natively, from the Mass Ledger's derived skin profile lambda(u) and
the corpus-licensed coupling class. ZERO adjustable parameters.

ORDER GATE (honesty): all numbers are computed and printed BEFORE any
comparison against the l* = 3 band. The comparison is the last section.

Substrate (all cited, none tuned):
  - Skin E-L first integral (Mass Ledger swing 4/22, booked functional):
        lambda(u) * (hbar^2/8m) * rho0 * u'^2 / u = rho0 * u * Delta(u)
    Delta(u) = tau_b u^{2/3} - (a_v+tau_b) u + a_v ;  lambda(u) = 1/(1+8u)
    a_v = 15.75 MeV, C = a_v + tau_b = 35.85 MeV (the crown identity)
  - hbar^2/m_N = 41.5 MeV fm^2 (named import, same as ledger)
  - r0 = 1.2 fm (named import, same as ledger)
  - Confined ladder (Inner Census, derived): rungs N = 2 n_r + l,
    per-class capacity (N+1)(N+2); spacing hw fixed by matching the filled
    ladder's <r^2> to the saturated core (3/5)R^2, R = r0 A^{1/3}.
  - Single-quantum rate well: V(r) = -C * u(r) — the gross contact column
    (discount proportional to meshed contact, capped at capacity), the
    two-column ledger's own E = V + T split. No depth was fitted.
  - Coupling class (corpus-licensed, GROUNDING sect.4): micro spin-orbit is
    the ONE (v/c)^2 time-cone anholonomy (Thomas core), strength set by the
    local rate gradient:
        DE(l) = (2l+1)/2 * (hbar c)^2/(2 (m c^2)^2) * <(1/r) dV/dr>
    evaluated on the intruder candidate orbital (n_r = 0, l), against the
    DERIVED skin gradient. Sign of the bare anholonomy: co RAISED
    (attractive well). The App-5 channel rule (ratified) says co DEEPER.
    Both booked; see report.
  - Comparison candidate (identification grade, NOT derived, printed for
    the record only): the lock-patch cone anholonomy credit h = 0.95 MeV
    (Mass Ledger swings 32-36) applied per winding step:
    kappa_E(l) = h (2l+1)/2.
"""
import math

# ---- constants (named ledger imports) ----
HBAR2_OVER_M = 41.5          # MeV fm^2
K8 = HBAR2_OVER_M / 8.0      # hbar^2/8m = 5.1875 MeV fm^2 (ledger quotes 5.184)
HBARC = 197.327              # MeV fm
MC2 = 938.92                 # MeV
R0 = 1.2                     # fm
AV = 15.75                   # MeV
CCAP = 35.85                 # MeV  (C = a_v + tau_b, crown identity)
TAUB = CCAP - AV             # 20.10 MeV
RHO0 = 3.0 / (4.0 * math.pi * R0**3)   # fm^-3, consistent with R = r0 A^(1/3)
H_ANHOL = 0.95               # MeV, lock-patch cone anholonomy credit (comparison only)

def lam(u): return 1.0 / (1.0 + 8.0 * u)

def Delta(u):
    return TAUB * u**(2.0/3.0) - (AV + TAUB) * u + AV

# ---- 1. skin profile u(x), x measured outward from the saturation kink ----
# dx/du = - sqrt(K8 * lam(u)) / (u * sqrt(Delta(u)))
def build_profile(n=40000, u_min=1e-6, u_max=1.0 - 1e-9):
    # integrate in s = log(u) for tail accuracy: dx/ds = -sqrt(K8*lam/Delta)
    s0, s1 = math.log(u_max), math.log(u_min)
    xs, us = [0.0], [u_max]
    s, x = s0, 0.0
    ds = (s1 - s0) / n
    def f(s):
        u = math.exp(s)
        d = Delta(u)
        if d <= 0: d = 1e-12
        return -math.sqrt(K8 * lam(u) / d)
    for i in range(n):
        k1 = f(s); k2 = f(s + ds/2); k3 = f(s + ds/2); k4 = f(s + ds)
        x += ds * (k1 + 2*k2 + 2*k3 + k4) / 6.0
        s += ds
        xs.append(x); us.append(math.exp(s))
    return xs, us

XS, US = build_profile()

def u_of_x(x):
    if x <= 0: return 1.0
    if x >= XS[-1]: return US[-1] * math.exp(-(x - XS[-1]) / math.sqrt(K8 / AV))
    lo, hi = 0, len(XS) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if XS[mid] < x: lo = mid
        else: hi = mid
    t = (x - XS[lo]) / (XS[hi] - XS[lo])
    return US[lo] * (1 - t) + US[hi] * t

ELL_TAIL = math.sqrt(K8 / AV)
# skin surface-width check: 10-90 fall
def x_at_u(target):
    lo, hi = 0.0, XS[-1]
    for _ in range(80):
        mid = (lo + hi) / 2
        if u_of_x(mid) > target: lo = mid
        else: hi = mid
    return (lo + hi) / 2

# ---- 2. place the profile: core edge Rc such that total census = A ----
def census_integral(Rc, dx=0.005, xmax=None):
    if xmax is None: xmax = XS[-1] + 8 * ELL_TAIL
    s = Rc**3 / 3.0
    x = 0.0
    while x < xmax:
        u1 = u_of_x(x); u2 = u_of_x(x + dx)
        r1 = Rc + x; r2 = Rc + x + dx
        s += 0.5 * (u1 * r1 * r1 + u2 * r2 * r2) * dx
        x += dx
    return 4 * math.pi * RHO0 * s

def solve_Rc(A):
    lo, hi = 0.1, R0 * A**(1/3.0) + 3.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if census_integral(mid) < A: lo = mid
        else: hi = mid
    return (lo + hi) / 2

# ---- 3. confined ladder: exact <r^2> match -> b^2, hw ----
def rung_cap(N): return (N + 1) * (N + 2)

def ladder_msr_coeff(nq):
    """Fill nq quanta (one class) up the ladder; return sum of (N+3/2) weights."""
    tot, N, w = 0, 0, 0.0
    while tot < nq:
        take = min(rung_cap(N), nq - tot)
        w += take * (N + 1.5)
        tot += take
        N += 1
    return w  # <r^2> = b^2 * w / nq

def hw_of_A(A):
    nq = A // 2  # per class, symmetric two-class bookkeeping
    R = R0 * A**(1/3.0)
    msr = 0.6 * R * R                      # (3/5) R^2
    b2 = msr * nq / ladder_msr_coeff(nq)   # fm^2
    return HBAR2_OVER_M / b2, math.sqrt(b2)  # hw (MeV), b (fm)

# ---- 4. orbital average of (1/r) dV/dr on (n_r=0, l) ----
def orbital_avg(l, b, Rc, dr=0.005):
    """<(1/r) dV/dr> with V(r) = -C u(r - Rc); density ~ r^{2l+2} exp(-r^2/b^2)."""
    num = den = 0.0
    rmax = Rc + XS[-1] + 8 * ELL_TAIL
    r = dr
    h = 1e-4
    while r < rmax:
        w = math.exp((2 * l + 2) * math.log(r) - r * r / (b * b))
        x = r - Rc
        dVdr = -CCAP * (u_of_x(x + h) - u_of_x(x - h)) / (2 * h)  # = +C|u'| in skin
        num += w * dVdr / r * dr
        den += w * dr
        r += dr
    return num / den

# ---- 5. the run ----
def cum(N):
    return sum(rung_cap(k) for k in range(N + 1))

print("=" * 78)
print("EYEWALL PRICING — BLIND RUN (numbers first, comparison last)")
print("=" * 78)
print(f"substrate: a_v={AV}, C={CCAP}, tau_b={TAUB:.2f} MeV; hbar^2/8m={K8:.4f} MeV fm^2")
print(f"skin profile from first integral: tail ell = {ELL_TAIL:.4f} fm "
      f"(booked 0.574; measured diffuseness 0.55±0.06)")
x10, x90 = x_at_u(0.90), x_at_u(0.10)
print(f"skin 90->10 width = {x90 - x10:.3f} fm; full profile extent {XS[-1]:.2f} fm to u=1e-6")
print()
print("Per intruder winding l = N (candidate orbital n_r=0, l), priced at the")
print("two bracketing symmetric censuses: A_lo = 2*cum(N-1) (rung N opens),")
print("A_hi = 2*cum(N) (rung N complete). kappa = FULL co/counter split.")
print()
hdr = (f"{'l':>2} {'A':>4} {'hw=Dnu':>7} {'b':>5} {'Rc':>5} {'<(1/r)dV/dr>':>13} "
       f"{'kappa_A':>8} {'k_A/Dnu':>8} {'kappa_E':>8} {'k_E/Dnu':>8}")
print(hdr)
results = {}
for l in range(1, 9):
    N = l
    for tag, A in (("lo", 2 * cum(N - 1) if N >= 1 else 4), ("hi", 2 * cum(N))):
        hw, b = hw_of_A(A)
        Rc = solve_Rc(A)
        avg = orbital_avg(l, b, Rc)
        kA = (2 * l + 1) / 2.0 * HBARC**2 / (2 * MC2**2) * avg   # (v/c)^2 anholonomy class
        kE = H_ANHOL * (2 * l + 1) / 2.0                          # comparison candidate
        results[(l, tag)] = (A, hw, kA / hw, kE / hw)
        print(f"{l:>2} {A:>4} {hw:>7.3f} {b:>5.2f} {Rc:>5.2f} {avg:>13.4f} "
              f"{kA:>8.3f} {kA/hw:>8.4f} {kE:>8.3f} {kE/hw:>8.4f}")
print()

# threshold under each pricing: smallest l with kappa/Dnu >= 1 (take A_lo..A_hi band)
def threshold(idx):
    ls = []
    for tag in ("lo", "hi"):
        t = None
        for l in range(1, 200):
            if l <= 8:
                r = results[(l, tag)][idx]
            else:
                # extrapolate: ratio ~ const * (2l+1) using l=8 value
                r = results[(8, tag)][idx] * (2 * l + 1) / 17.0
            if r >= 1.0:
                t = l
                break
        ls.append(t)
    return ls

tA = threshold(2)
tE = threshold(3)
print(f"BLIND RESULT — licensed (v/c)^2 anholonomy pricing: kappa/Dnu at l=3 = "
      f"{results[(3,'lo')][2]:.4f}..{results[(3,'hi')][2]:.4f}")
print(f"  threshold l* (first kappa/Dnu >= 1): {tA[0]}..{tA[1]}")
print(f"  bare anholonomy sign: co-turning RAISED (attractive well) — "
      f"App-5 channel rule says co DEEPER; bare sign anti-record.")
print(f"COMPARISON candidate (lock-patch cone credit, NOT derived): kappa/Dnu at l=3 = "
      f"{results[(3,'lo')][3]:.4f}..{results[(3,'hi')][3]:.4f}; threshold l* = {tE[0]}..{tE[1]}")
print()

# what the RECORD requires of any pricing (for the report, not a target of the run):
for tag in ("lo",):
    A, hw, _, _ = results[(3, tag)]
    print(f"record-implied requirement at the 28 seam: split ~ Dnu = {hw:.2f} MeV at l=3, A~{A}")

print()
print("=" * 78)
print("COMPARISON (the gate, applied last)")
print("=" * 78)
okA = (tA[0] == 3 and tA[1] == 3)
okE = (tE[0] == 3 and tE[1] == 3)
print(f"licensed pricing selects l* = 3 exactly: {'PASS' if okA else 'FAIL'}")
print(f"comparison candidate selects l* = 3 exactly: {'PASS' if okE else 'FAIL'}")
if not okA:
    m = 1.0 / results[(3, 'lo')][2]
    print(f"licensed pricing misses the l*=3 band by a factor ~{m:.0f} at l=3 "
          f"(and carries the anti-record bare sign).")
print()
print("VERDICT: see report. Growth ~ (2l+1) against constant spacing (structural")
print("threshold existence) is confirmed by both pricings; the LOCATION under every")
print("corpus-licensed zero-knob pricing is far above l*=3. The Inner Census")
print("falsifier condition ('threshold at l* != 3') is TRIGGERED for the licensed")
print("coupling class.")

# ============================================================================
# PHASE-3 CANDIDATE (added 2026-09-06, ratified): the CAPPED TRANSPORT LAW.
# The gate arithmetic above is FROZEN and untouched; this section only adds a
# new candidate and prices it through the same criteria (threshold = first l
# with kappa/Dnu >= 1; unit descent = ratio < 2 for realized l = 3..6; band
# l* = 3). The licensed-class FAIL above stands as the historical control.
#
# Derivation chain (identification grade — three joints NAMED, not forced):
#  (1) Transactions bank per completed closure (whole-turn primacy); a
#      winding-l mode completes l phase closures per orbit.
#  (2) Each closure is co-signed through ONE lock; the co-signing lock's
#      patch cone is Omega_w = 4pi/z_c EXACTLY (Mass Ledger swing 23).
#      [JOINT J1: one closure wraps one patch cone]
#  (3) Sense-pair transport around the cone splits the pair's phase by
#      Omega_w  =>  per-step price (2/z_c) * Dnu.
#  (4) INDEPENDENCE CAP (the ledger's own cycle-space theorem, swings 28-31):
#      the cones tile the quantum's sphere exactly (z_c * Omega_w = 4pi);
#      wrapping all of them is the trivial cycle, so at most (z_c - 1)
#      independent wraps bank per orbit: cap = 2(1 - 1/z_c) < 2 for ANY z_c
#      — unit descent is structural.  [JOINT J2: intrinsic tiling; the
#      integer-seam variant 8/z_c = 1.62 noted, gates insensitive]
#  (5) Sign: co-turning sense BANKS the transport credit (deeper); counter
#      pays. Agrees with the App-5 channel rule and the record.
#      [JOINT J3: bank<->co orientation]
#
#      LAW:  kappa(l)/Dnu = min( 2l/z_c , 2(1 - 1/z_c) )
# ============================================================================

print()
print("=" * 78)
print("CAPPED TRANSPORT LAW (ratified 2026-09-06; identification grade, joints J1-J3)")
print("=" * 78)
ZC_BAND2 = [("tangent-cone chain", 4.940), ("realized packing", 4.846),
            ("energy route", 4.780), ("band top", 5.000)]
REALIZED = [3, 4, 5, 6]
all_pass = True
for zname, zc in ZC_BAND2:
    law = lambda l: min(2.0 * l / zc, 2.0 * (1.0 - 1.0 / zc))
    row = "  ".join(f"l={l}:{law(l):.3f}" for l in range(1, 9))
    lstar = next(l for l in range(1, 60) if law(l) >= 1.0)
    unit = all(law(l) < 2.0 for l in REALIZED)
    persist = all(law(l) >= 1.0 for l in range(3, 20))
    r184 = law(7)
    ok = (lstar == 3) and unit and persist and (1.0 <= r184 < 2.0)
    all_pass &= ok
    print(f"z_c = {zc} ({zname}):")
    print(f"  {row}")
    print(f"  l* = {lstar} | unit descent l=3..6: {'OK' if unit else 'BROKEN'} | "
          f"persistence: {'OK' if persist else 'NO'} | 184 (l=7): {r184:.3f} -> "
          f"{'UNIT descent' if 1.0 <= r184 < 2.0 else 'other'} | "
          f"gate: {'PASS' if ok else 'FAIL'}")
print()
print(f"structural window for l* = 3: 4 < z_c <= 6 (derived band [4.78, 5.00] strictly inside)")
print(f"unit descent: THEOREM (cap 2(1-1/z_c) < 2 for any finite z_c, all l)")
print(f"persistence:  THEOREM (cap >= 1 iff z_c >= 2)")
print(f"CAPPED LAW across the derived z_c band: {'PASS on all gates' if all_pass else 'FAIL'}")
print()
print("Booked alongside (standalone zero-knob exclusion): the composed-centrifugal")
print("backreaction (j replacing l in the radial books) gives kappa/Dnu =")
print("(2l+1)/(2l+3) < 1 for ALL l — sense composition alone can never re-file a")
print("multiplet. Tagged extant tension, falsifiable: the law's split (~1.2 hw at")
print("l=3) is ~2x extant-extracted spectroscopic 1f splittings; this page prices")
print("seam bookkeeping, not level spectroscopy.")

# ============================================================================
# PHASE-5 CONTROL (added 2026-09-06, ratified): SMOOTH (CONTINUUM) TRANSPORT.
# Gate arithmetic untouched. This is the law's SECOND control, beside the
# licensed (v/c)^2 miss above.
#
# The orbit plane passes through the core, so the smoothly-transported
# co-signing direction traces a GREAT CIRCLE on the direction sphere per
# circuit: enclosed solid angle = 2pi exactly — every orbit, every l, every
# eccentricity (pericenter closure never enters). Spin-1/2 pair split under
# smooth transport = hbar * 2pi / T_orb = hw: EXACTLY 1.00 Dnu at every
# l >= 1, l-INDEPENDENT.
# ============================================================================

print()
print("=" * 78)
print("SMOOTH-TRANSPORT CONTROL (continuum reading; ratified 2026-09-06)")
print("=" * 78)
smooth = lambda l: 1.00  # 2pi great-circle sweep -> split = hw at every l >= 1
row = "  ".join(f"l={l}:{smooth(l):.2f}" for l in range(1, 8))
print(f"  kappa/Dnu: {row}")
lstar_s = next(l for l in range(1, 60) if smooth(l) >= 1.0)
print(f"  threshold l* = {lstar_s} -> EVERY rung's co-multiplet intrudes from l = 1:")
print(f"  the 8 and 20 seams would not exist. THE RECORD'S INTACT 8/20 SEAMS")
print(f"  POSITIVELY EXCLUDE CONTINUUM TRANSPORT: FAIL by record.")
print(f"  Discrete per-closure accounting is DATA-SELECTED, not assumed —")
print(f"  'integrate, don't differentiate' at the transport level.")
print()
print("J1 status after the Phase-5 theorem work (see inner-census.html §6b):")
print("  J1a PROVEN (kinematic: l whole phase turns per circuit; n_r=0 modes —")
print("      every realized intruder — carry no other closure entries);")
print("  J1b RATIFIED ONTOLOGY (one entry per turn, one co-signature per entry,")
print("      footprint Omega_w = 4pi/z_c);")
print("  J1c ONE NAMED IDENTIFICATION (wrap-per-entry: integer wraps forced by")
print("      whole-turn discipline; cluster-style fractions are sum-preserving")
print("      latency, cancelling in the orbit sum).")
print("  Fraction-robustness window the gates themselves measure: f in")
print("  (z_c/6, z_c/4] = (0.823, 1.235]; theorem value f = 1 central; no z_c")
print("  narrowing results (non-result booked).")
print("  GRADE: derivation with one named identification (was: identification,")
print("  three named joints).")
