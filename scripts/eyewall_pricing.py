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
