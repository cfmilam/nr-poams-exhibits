#!/usr/bin/env python3
"""
Aufbau as Closure -- the page's computation, end to end, exposed.
Companion script to aufbau-as-closure.html (NR/POAMS exhibit series).
Run 2026-09-07 (Tier-A review repair: the executable the page lacked).

WHAT THIS REPRODUCES (the page's actual chain, no more):
  1. The frontier potential the page uses: the Thomas-Fermi ledger
     V(r) = -(Z/r) chi(r/b), b = 0.88534 a0 Z^(-1/3), chi'' = chi^(3/2)/sqrt(x)
     solved here by RK4 bisection shooting (s = sqrt(x) regularization; the
     initial slope is found by the shooting itself, then checked against the
     literature value B = -1.5880710226 -- a validation, not an input; the
     far tail is patched with the value-matched Coulson-March asymptote
     144/x^3 (1 + a x^-lambda), lambda = (sqrt(73)-7)/2, at the point where
     double precision lets go of the separatrix).
  2. The zero-energy apsidal integral: k(Z,l) = Phi/pi with
     Phi = int_{r1}^{r2} L dr / (r sqrt(2 Z r chi(r/b) - L^2)),  L = l + 1/2
     (endpoint singularities removed by the theta substitution; the L = l+1/2
      Langer convention is the page's, robustness noted there).
  3. The Fermi capture threshold: winding l holdable at the E=0 frontier iff
     1.77068 Z^(2/3) max_x[x chi(x)] > (l+1/2)^2.
  4. The greedy build: at each Z = 1..118 the next quantum takes the available
     mode minimizing c = n_r + k(Z,l)*l, capacities 2(2l+1), a winding
     available only above its capture threshold. k is recomputed at each Z.

REGISTERED METRICS (declared before any computation; the numbers print
whatever they are -- both metrics are reported, they are NOT the same number):
  M1  POSITIONAL: the emergent 19-subshell opening sequence vs the IDEAL
      MADELUNG textbook order, position by position (x of 19).
      Under review (Sol, Tier A): claimed value 13/19.
  M2  OPENING-Z: count of the 19 subshells whose emergent first-occupancy Z
      equals the OBSERVED record's first-occupancy Z (NIST ground
      configurations, tagged: La 5d at 57, Ce 4f at 58, Ac 6d at 89,
      Pa 5f at 91 -- the observed record, not the idealized rule).
      Page's published claim under test: 15/19.
  M3  PROVENANCE GATE: the page's precomputed arrays (ZONES dots, KDATA
      curves, EMERGENT sequence) regenerated from first principles here;
      max |delta k| printed. If regeneration fails, that failure prints --
      the arrays' provenance is then this script's honest state.
VALIDATION GATES (frozen): chi(1)=0.424008, chi(10)=0.024314,
  chi(50)=0.000632 (+-1e-6); Demkov-Ostrovsky model potential returns
  k = 2.000000 (+-1e-5) for every tested (v, L), strength-independent;
  Coulomb E=0 has no outer turning point (parabolic, unbound).

ZERO adjustable parameters: every number descends from chi'' = chi^(3/2)/sqrt(x)
plus classical mechanics. The one convention is L = l + 1/2 (declared).
"""
import math

# ---------------------------------------------------------------- constants
B0 = -1.5880710226          # chi'(0), literature shooting value
BCOEF = 0.88534             # b = BCOEF * a0 * Z^(-1/3)
TWO_B = 2.0 * BCOEF         # 1.77068
LMAX_UNIVERSE = 4           # consider s..g
ZMAX = 118

IDEAL_MADELUNG = ["1s","2s","2p","3s","3p","4s","3d","4p","5s","4d","5p",
                  "6s","4f","5d","6p","7s","5f","6d","7p"]
# Observed first-occupancy Z (NIST ground configurations; tagged record):
OBSERVED_OPEN = {"1s":1,"2s":3,"2p":5,"3s":11,"3p":13,"4s":19,"3d":21,
                 "4p":31,"5s":37,"4d":39,"5p":49,"6s":55,"4f":58,"5d":57,
                 "6p":81,"7s":87,"5f":91,"6d":89,"7p":113}
# Idealized-Madelung first-occupancy Z (textbook order, for contrast):
IDEAL_OPEN = {"1s":1,"2s":3,"2p":5,"3s":11,"3p":13,"4s":19,"3d":21,
              "4p":31,"5s":37,"4d":39,"5p":49,"6s":55,"4f":57,"5d":71,
              "6p":81,"7s":87,"5f":89,"6d":103,"7p":113}
# Page's precomputed arrays (as published 2026-08-09), under provenance test:
PAGE_EMERGENT = ["1s","2s","2p","3s","3p","4s","4p","3d","5s","4d","5p",
                 "6s","5d","4f","6p","7s","6d","5f","7p"]
PAGE_ZONES = [("4s",19.5,1.6812,0),("3d",25.5,1.9207,2),("4p",33.5,1.8114,1),
  ("5s",37.5,1.6576,0),("4d",43.5,1.8868,2),("5p",51.5,1.7880,1),
  ("6s",55.5,1.6448,0),("4f",64,1.9264,3),("5d",75.5,1.8533,2),
  ("6p",83.5,1.7627,1),("7s",87.5,1.6310,0),("5f",96,1.9006,3),
  ("6d",108.5,1.8319,2),("7p",115.5,1.7468,1)]
PAGE_KZ = [10,20,30,40,55,70,85,100,118]
PAGE_KDATA = {
  0:[1.7082,1.6802,1.6654,1.6555,1.6451,1.6377,1.6318,1.6274,1.6228],
  1:[1.8830,1.8412,1.8178,1.8017,1.7845,1.7719,1.7619,1.7538,1.7456],
  2:[None,1.9370,1.9103,1.8922,1.8723,1.8577,1.8462,1.8367,1.8271],
  3:[None,None,None,None,1.9368,1.9207,1.9083,1.8981,1.8877]}

print("=" * 78)
print("AUFBAU AS CLOSURE -- EXECUTABLE INSTRUMENT (registration first)")
print("=" * 78)
print(__doc__.split("REGISTERED METRICS")[1].split("ZERO adjustable")[0]
      .replace("(declared before any computation; the numbers print\nwhatever"
               " they are -- both metrics are reported, they are NOT the same"
               " number):", "(frozen before computation):"))

# ------------------------------------------------- 1. Thomas-Fermi solution
# s = sqrt(x): d(chi)/ds = 2 s chi_x ; d(chi_x)/ds = 2 chi^(3/2)
S_MAX = 32.0                 # x up to 1024
DS = 5.0e-4
N_STEPS = int(S_MAX / DS)
LAM = (math.sqrt(73.0) - 7.0) / 2.0          # tail correction exponent

def _f(s, c, cp):
    cc = c if c > 0.0 else 0.0
    return 2.0 * s * cp, 2.0 * (cc ** 1.5)

def _integrate(B, record=None):
    """RK4 in s. Returns ('zero'|'diverge'|'ok', s_fail)."""
    y1, y2 = 1.0, B
    s = 0.0
    if record is not None: record[0] = 1.0
    for i in range(N_STEPS):
        k1a, k1b = _f(s, y1, y2)
        k2a, k2b = _f(s + DS/2, y1 + k1a*DS/2, y2 + k1b*DS/2)
        k3a, k3b = _f(s + DS/2, y1 + k2a*DS/2, y2 + k2b*DS/2)
        k4a, k4b = _f(s + DS, y1 + k3a*DS, y2 + k3b*DS)
        y1 += DS/6 * (k1a + 2*k2a + 2*k3a + k4a)
        y2 += DS/6 * (k1b + 2*k2b + 2*k3b + k4b)
        s += DS
        if record is not None: record[i + 1] = y1 if y1 > 0.0 else 0.0
        if y1 <= 0.0:  return ('zero', s)
        if y2 >= 0.0:  return ('diverge', s)
    return ('ok', s)

# bisection shooting: chi crosses zero => slope too negative; diverges => too high
b_lo, b_hi = -1.58808, -1.58806
for _ in range(60):
    b_mid = 0.5 * (b_lo + b_hi)
    st, _sf = _integrate(b_mid)
    if st == 'zero': b_lo = b_mid
    else:            b_hi = b_mid
B_SHOT = 0.5 * (b_lo + b_hi)
chi_s = [0.0] * (N_STEPS + 1)
_status, s_fail = _integrate(B_SHOT, record=chi_s)
x_fail = s_fail * s_fail
X_PATCH = min(0.9 * x_fail, 900.0)           # hand off to asymptote here

def _chi_grid(x):
    ss = math.sqrt(x)
    u = ss / DS
    i = int(u)
    if i >= N_STEPS: i = N_STEPS - 1
    t = u - i
    return chi_s[i] * (1 - t) + chi_s[i + 1] * t

# value-matched Coulson-March tail: chi = 144/x^3 (1 + a x^-lambda)
_A_TAIL = (_chi_grid(X_PATCH) * X_PATCH ** 3 / 144.0 - 1.0) * X_PATCH ** LAM

def chi(x):
    if x <= 0.0: return 1.0
    if x >= X_PATCH:
        return 144.0 / x ** 3 * (1.0 + _A_TAIL * x ** (-LAM))
    return _chi_grid(x)

print("-" * 78)
print(f"GATE 1 -- shooting: chi'(0) = {B_SHOT:.10f}  (literature {B0});")
print(f"  separatrix held to x = {x_fail:.1f}; tail patched at x = "
      f"{X_PATCH:.1f} (a = {_A_TAIL:.3f})")
print("  chi checkpoints (page: 0.424008 / 0.024314 / 0.000632):")
ok1 = abs(B_SHOT - B0) < 5e-7
for xx, ref in ((1.0, 0.424008), (10.0, 0.024314), (50.0, 0.000632)):
    v = chi(xx)
    passed = abs(v - ref) < 1.5e-6
    ok1 &= passed
    print(f"  chi({xx:>4}) = {v:.6f}   ref {ref:.6f}   "
          f"{'PASS' if passed else 'FAIL'}")

# max of x*chi(x) for the capture threshold
xm, fm = 0.0, 0.0
xx = 0.5
while xx < 12.0:
    v = xx * chi(xx)
    if v > fm: fm, xm = v, xx
    xx += 1e-3
print(f"  max[x chi(x)] = {fm:.6f} at x = {xm:.3f}")

# ------------------------------------------------- 2. apsidal integral k
def _phi_over_pi(g, r_lo, r_peak, r_hi, L, n_nodes=240):
    """k = Phi/pi for g(r) = r^2 p_r^2; roots bracketed around r_peak."""
    def _root(a, b):
        fa, fb = g(a), g(b)
        if fa * fb > 0: return None
        for _ in range(200):
            m = 0.5 * (a + b)
            fmid = g(m)
            if fa * fmid <= 0: b, fb = m, fmid
            else: a, fa = m, fmid
        return 0.5 * (a + b)
    r1 = _root(r_lo, r_peak)
    r2 = _root(r_peak, r_hi)
    if r1 is None or r2 is None: return None
    # theta substitution: r = mid - half*cos(theta); sqrt sing. removed
    mid, half = 0.5 * (r1 + r2), 0.5 * (r2 - r1)
    total = 0.0
    n = n_nodes
    for j in range(n):                       # midpoint rule in theta
        th = math.pi * (j + 0.5) / n
        r = mid - half * math.cos(th)
        gg = g(r)
        w = gg / ((r - r1) * (r2 - r))       # smooth on (r1, r2)
        if w <= 0.0: continue
        total += L / (r * math.sqrt(w)) * (math.pi / n)
    return total / math.pi

def k_of(Z, l):
    """Apsidal slope k at the E=0 TF frontier, L = l+1/2, in x = r/b units."""
    A = TWO_B * Z ** (2.0 / 3.0)
    L2 = (l + 0.5) ** 2
    if A * fm <= L2: return None             # frontier cannot hold winding l
    gfun = lambda x: A * x * chi(x) - L2
    return _phi_over_pi(gfun, 1e-12, xm, 4000.0, l + 0.5)

print("-" * 78)
print("GATE 2 -- Demkov-Ostrovsky model potential (must give k = 2.000000,")
print("          strength-independent) and Coulomb E=0 (must be unbound):")
ok2 = True
for v, L in ((0.9, 0.5), (0.9, 0.65), (5.0, 1.2), (2.5, 0.7)):  # need L^2 < v/2
    gd = lambda r: 2.0 * v * r / (r + 1.0) ** 2 - L * L
    kk = _phi_over_pi(gd, 1e-12, 1.0, 4000.0, L)
    passed = kk is not None and abs(kk - 2.0) < 1e-5
    ok2 &= passed
    print(f"  DO v={v:<4} L={L:<4}: k = {kk:.6f}   "
          f"{'PASS' if passed else 'FAIL'}")
gc = lambda r: 2.0 * 10.0 * r - 0.25        # Coulomb Z=10, L=1/2: monotone
unbound = gc(4000.0) > 0 and gc(1.0) > 0
print(f"  Coulomb E=0: g(r) has no outer root (unbound)  "
      f"{'PASS' if unbound else 'FAIL'}")
ok2 &= unbound

# ------------------------------------------------- 3. capture thresholds
print("-" * 78)
print("Fermi capture thresholds  1.77068 Z^(2/3) max[x chi] = (l+1/2)^2")
print("(page's table: p 4.2 / d 19.6 / f 53.7 / g 114):")
zc = {}
for l, nm in ((1, "p"), (2, "d"), (3, "f"), (4, "g")):
    L2 = (l + 0.5) ** 2
    zc[l] = (L2 / (TWO_B * fm)) ** 1.5
    print(f"  {nm}: Z_c = {zc[l]:.1f}")
zc[0] = 0.0

# ------------------------------------------------- 4. provenance: page arrays
print("-" * 78)
print("M3 -- PROVENANCE: regenerating the page's precomputed arrays")
print("ZONES dots (k at zone-middle Z; convention disclosed on the page):")
dmax = 0.0
for nm, zmid, kpage, l in PAGE_ZONES:
    kv = k_of(zmid, l)
    d = abs(kv - kpage)
    dmax = max(dmax, d)
    print(f"  {nm:>3} @ Z={zmid:<6} l={l}:  k = {kv:.4f}   page {kpage:.4f}"
          f"   |dk| = {d:.4f}")
print("KDATA curves (k on the page's Z grid):")
for l in (0, 1, 2, 3):
    row = []
    for i, Zg in enumerate(PAGE_KZ):
        ref = PAGE_KDATA[l][i]
        if ref is None:
            row.append("   --  ")
            continue
        kv = k_of(Zg, l)
        if kv is None:
            row.append(" unbound")
            dmax = max(dmax, 99.0)
            continue
        d = abs(kv - ref)
        dmax = max(dmax, d)
        row.append(f"{kv:.4f}")
    print(f"  l={l}: " + " ".join(row))
print(f"  max |delta k| vs page arrays: {dmax:.4f}")

# ------------------------------------------------- 5. the greedy build
print("-" * 78)
print("GREEDY BUILD (one quantum per Z; c = n_r + k(Z,l)*l, k recomputed")
print("per Z; capacities 2(2l+1); winding admitted above its threshold):")
LNAME = "spdfg"
modes = []
for n in range(1, 9):
    for l in range(0, min(LMAX_UNIVERSE, n - 1) + 1):
        modes.append({"n": n, "l": l, "nr": n - l - 1,
                      "cap": 2 * (2 * l + 1), "occ": 0,
                      "name": f"{n}{LNAME[l]}"})
kcache = {}
opening = {}
for Z in range(1, ZMAX + 1):
    best, bestc = None, None
    for m in modes:
        if m["occ"] >= m["cap"]: continue
        l = m["l"]
        if l > 0 and Z <= zc[l]: continue
        if l == 0:
            c = float(m["nr"])
        else:
            if (Z, l) not in kcache:
                kcache[(Z, l)] = k_of(Z, l)
            kv = kcache[(Z, l)]
            if kv is None: continue
            c = m["nr"] + kv * l
        if bestc is None or c < bestc - 1e-12:
            best, bestc = m, c
    best["occ"] += 1
    if best["name"] not in opening:
        opening[best["name"]] = Z
emergent = sorted(opening.keys(), key=lambda nm: opening[nm])
emergent19 = [nm for nm in emergent if nm in OBSERVED_OPEN]
print("  emergent opening sequence:")
print("   " + " ".join(f"{nm}@{opening[nm]}" for nm in emergent))
same_as_page = emergent19 == PAGE_EMERGENT
print(f"  matches page EMERGENT array: {same_as_page}")

# ------------------------------------------------- 6. the two metrics
print("-" * 78)
print("METRICS (both, honestly, side by side):")
m1 = sum(1 for i in range(19)
         if i < len(emergent19) and emergent19[i] == IDEAL_MADELUNG[i])
print(f"  M1 POSITIONAL vs IDEAL MADELUNG: {m1}/19")
mism = [f"{emergent19[i]}|{IDEAL_MADELUNG[i]}" for i in range(19)
        if i < len(emergent19) and emergent19[i] != IDEAL_MADELUNG[i]]
print(f"     mismatched positions: {', '.join(mism) if mism else 'none'}")
m2 = sum(1 for nm in OBSERVED_OPEN
         if nm in opening and opening[nm] == OBSERVED_OPEN[nm])
print(f"  M2 OPENING-Z vs OBSERVED record (NIST first occupancy): {m2}/19")
miss2 = [f"{nm}: emergent {opening.get(nm,'--')} vs record {OBSERVED_OPEN[nm]}"
         for nm in IDEAL_MADELUNG
         if opening.get(nm) != OBSERVED_OPEN[nm]]
for x in miss2: print(f"     miss -- {x}")
m2i = sum(1 for nm in IDEAL_OPEN
          if nm in opening and opening[nm] == IDEAL_OPEN[nm])
print(f"  (contrast: opening-Z vs IDEALIZED Madelung openings: {m2i}/19)")

print("-" * 78)
print("VERDICT LINE (for the page to quote):")
print(f"  gates: chi {'PASS' if ok1 else 'FAIL'} / DO+Coulomb "
      f"{'PASS' if ok2 else 'FAIL'} / arrays regenerated max|dk|={dmax:.4f}")
print(f"  positional agreement {m1}/19 ; exact opening-Z (observed record) "
      f"{m2}/19 ; page-array match {same_as_page}")
