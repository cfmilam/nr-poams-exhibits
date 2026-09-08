#!/usr/bin/env python3
"""
The Inner Census — counting verification.

Companion script to inner-census.html (NR/POAMS exhibit series).

Verifies two distinct arithmetic statements without promoting either one's
physical premises: (i) a confined harmonic ladder has cumulative closures
2, 8, 20 exactly; and (ii) the historical intruder bookkeeping reproduces
28, 50, 82, 126 (and 184) when its threshold, one-rung descent, and sign are
read from the measured record. Both native pricings of that intruder rule
were later falsified. The script remains a counting certificate, not a
derivation of the physical nuclear sequence.

The rule
--------
1. THE CONFINED LADDER (exact conditional arithmetic).
   Inside a specified saturated, flat-rate candidate core, the relevant closure-compatible
   geometry is the confined one: pericenter every half turn, so one
   radial coil trades for exactly two winding steps.  Rungs are
   labelled N = 2*n_r + l.  Rung N holds windings l = N, N-2, ..., 1|0,
   each with (2l+1) orientations, each orientation holding the binary
   pair of the quantum's own turn (the sense pair).  Per class:

       cap(N) = 2 * sum_{l = N, N-2, ...} (2l+1) = (N+1)(N+2)

   Cumulative: 2, 8, 20, 40, 70, 112, 168, ...

2. THE CO-TURN SPLIT (derived counting; sign channel-motivated).
   Composing the winding ladder (2l+1 orientations) with the binary
   sense pair yields exactly two combined ladders: the co-turning
   multiplet of size 2l+2 and the counter-turning multiplet of
   size 2l.   (2l+2) + (2l) = 2(2l+1).  Pure counting.

3. THE INTRUDER RULE (record-description; mechanism open).
   From the measured record the threshold is l* = 3: for every rung N >= 3 the
   top co-turning multiplet (l = N, size 2N+2) descends exactly one
   rung, landing just ABOVE the band below.  A closure gap always
   survives above the descended multiplet.  Below it: if the band
   underneath is a complete, closed rung (below threshold), the seam
   under the intruder survives too (the 20/28 double gap); if the
   band underneath is a broken remainder, the intruder joins its top
   and the vacated ladder seam is destroyed as a closure (40, 70,
   112 survive at most as the record's semi-magic seams).

Expected arithmetic: 2, 8, 20, 28, 50, 82, 126 (and conditional 184).
"""

L_STAR = 3  # the one imported integer: intrusion threshold (named, not derived)
N_MAX = 7   # build rungs 0..7 (rung 7 supplies the intruder for the 184 group)

TARGET = [2, 8, 20, 28, 50, 82, 126]


def rung_windings(N):
    """Windings on rung N of the confined ladder: l = N, N-2, ..., 1 or 0."""
    return list(range(N, -1, -2))


def rung_capacity(N):
    """Per-class capacity of rung N: 2 * sum(2l+1) over the rung's windings."""
    return 2 * sum(2 * l + 1 for l in rung_windings(N))


def multiplets(N):
    """All co/counter multiplets of rung N as (l, kind, size)."""
    out = []
    for l in rung_windings(N):
        out.append((l, "co", 2 * l + 2))       # co-turning: one step longer
        if l > 0:
            out.append((l, "counter", 2 * l))  # counter-turning
    return out


def check(label, got, want):
    ok = got == want
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {got}" + ("" if ok else f"  (expected {want})"))
    if not ok:
        raise SystemExit(f"VERIFICATION FAILED at: {label}")


print("== 1. The confined ladder (exact conditional arithmetic) ==")
caps = [rung_capacity(N) for N in range(N_MAX + 1)]
check("rung capacities (N+1)(N+2)", caps[:6], [2, 6, 12, 20, 30, 42])
check("closed form matches explicit count",
      caps, [(N + 1) * (N + 2) for N in range(N_MAX + 1)])

cum, tot = [], 0
for c in caps:
    tot += c
    cum.append(tot)
check("pure-ladder cumulative closures", cum[:6], [2, 8, 20, 40, 70, 112])
check("first three confined-ladder closures (given the premise)", cum[:3], TARGET[:3])

print("\n== 2. The co-turn split (derived counting) ==")
for N in range(N_MAX + 1):
    s = sum(size for _, _, size in multiplets(N))
    assert s == rung_capacity(N), f"multiplet sizes do not tile rung {N}"
print("  [PASS] co (2l+2) + counter (2l) multiplets exactly tile every rung 0..%d" % N_MAX)
check("top co-multiplet sizes 2N+2 (the intruders)",
      [2 * N + 2 for N in range(3, 8)], [8, 10, 12, 14, 16])

print("\n== 3. Imported intruder rule: assemble groups, read the gaps ==")
# groups = lists of multiplets between consecutive closure gaps,
# generated from the descent rule.  Start from the pure ladder
# (one group per rung), then move each rung-N co-multiplet (N >= l*)
# down one rung: it lands as its own group when the band below is a
# complete closed rung (N-1 < l*), else it joins the top of the
# group below (the remainder of rung N-1).
groups = [multiplets(N) for N in range(N_MAX + 1)]
descended = []
for N in range(L_STAR, N_MAX + 1):
    intruder = (N, "co", 2 * N + 2)
    groups[N].remove(intruder)
    descended.append((N, intruder))
for N, intruder in descended:
    if N - 1 < L_STAR:
        # band below is a complete closed rung: sit in the gap, own group
        groups.insert(N, [intruder])
    else:
        # band below is a remainder: join its top
        target = next(g for g in groups if (N - 1, "counter", 2 * (N - 1)) in g)
        target.append(intruder)

closures, tot = [], 0
for g in groups:
    tot += sum(size for _, _, size in g)
    closures.append(tot)
# The last group (rung N_MAX's remainder) is incomplete bookkeeping — its
# closing intruder would come from rung N_MAX+1 — drop the trailing entry.
closures = closures[:-1]

print("  groups (per class):")
names = "spdfghijk"
for i, g in enumerate(groups[:-1]):
    lab = " + ".join(f"{l}{names[l]}·{k}({s})" for l, k, s in g)
    print(f"    up to {closures[i]:>3}: {lab}")

check("closure sequence", closures[:7], TARGET)
check("conditional next closure (rung-6 remainder + imported rung-7 intruder)", closures[7], 184)

print("\n== 4. Closed form for the reordered numbers ==")
# magic(N) = cumulative ladder through rung N + intruder from rung N+1
#          = (N+1)(N+2)(N+3)/3 + (2N+4) = (N+2)(N^2+4N+9)/3,  N >= 2
cf = [(N + 2) * (N * N + 4 * N + 9) // 3 for N in range(2, 7)]
check("(N+2)(N²+4N+9)/3 for N=2..6", cf, [28, 50, 82, 126, 184])

print("\nALL CHECKS PASS — the stated counting rule reproduces 2, 8, 20, 28, 50, 82, 126 exactly")
print("GRADE: 2/8/20 are exact given the confined-ladder premise; 28+ are a")
print("record-description using imported threshold/descent/sign rules whose native pricings failed.")
