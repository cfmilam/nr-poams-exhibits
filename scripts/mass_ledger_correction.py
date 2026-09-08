#!/usr/bin/env python3
"""Mass Ledger finite-cluster correction — registered audit tier (2026-09-08).

This instrument does one bounded job: apply the intrinsic Robertson bound with
the finite-n factor squared, reproduce the corrected light-cluster lower-bound
quotients, and expose the blast radius of the former flat-contact inference.

The arithmetic is printed in the ledger's native whole-turn rate first.  MeV is
retained only in tagged parentheses so the result can be checked against the
record and the original swing-6 registration.
"""

from __future__ import annotations

MEV_TO_ZHZ = 0.2417989242
ZHZ_TO_MEV = 1.0 / MEV_TO_ZHZ
HBAR2_OVER_M = 41.47 * MEV_TO_ZHZ  # native ZHz fm^2 (record: 41.47 MeV fm^2)
A_C = 0.72 * MEV_TO_ZHZ            # native ZHz (record: 0.72 MeV)
GAMMA_CENTRAL = 0.85       # named channel input; not re-derived here
GAMMA_BAND = (0.80, 0.85, 0.90)

# label, Z, N, measured total binding (record MeV), inherited operational
# rms-radius proxy for sqrt(<x_i^2>) (fm), n, and the former channel-weight
# bookkeeping (strong, weak). The correction does not upgrade the radius proxy.
CLUSTERS = (
    ("d",     1, 1,  2.2245662, (1.950, 1.965, 1.980), 2, (1, 0)),
    ("He-3",  2, 1,  7.7180430, (1.750, 1.785, 1.820), 3, (1, 2)),
    ("t",     1, 2,  8.4817962, (1.540, 1.590, 1.680), 3, (1, 2)),
    ("alpha", 2, 2, 28.2956624, (1.450, 1.452, 1.480), 4, (2, 4)),
)

EXPECTED_ZHZ = {
    "d": 1.9987,
    "He-3": 2.5290,
    "t": 2.9631,
    "alpha": 3.5371,
}


def record_mev(native_zhz: float) -> float:
    """Tagged extant-record orientation for a native whole-turn rate."""
    return native_zhz * ZHZ_TO_MEV


def intrinsic_tax_lower_bound(n: int, rms_radius_fm: float) -> float:
    """Per-quantum Robertson lower bound, returned as native ZHz."""
    finite_n = 1.0 - 1.0 / n
    return (9.0 / 8.0) * HBAR2_OVER_M * finite_n**2 / rms_radius_fm**2


def strain(Z: int, A: int) -> float:
    return A_C * Z * (Z - 1) / A ** (1.0 / 3.0)


def lower_bound_quotient(cluster, radius: float, gamma: float) -> float:
    _, Z, N, binding, _, n, (strong, weak) = cluster
    A = Z + N
    binding_rate = binding * MEV_TO_ZHZ  # measured record imported, then cashed into rate
    gross_lower_bound = binding_rate + A * intrinsic_tax_lower_bound(n, radius) + strain(Z, A)
    return gross_lower_bound / (strong + weak * gamma)


print("MASS LEDGER CORRECTION — registration before confrontation")
print("  intrinsic variables: x_i = r_i - R_cm; q_i = p_i - P/n")
print("  [x_iα,q_iα] = i*hbar*(1-1/n)")
print("  registered correction: Robertson squares that factor")
print("  no radius, gamma, or channel weight is retuned")
print()
print("nu_tau,min = (9/8h)(hbar^2/m)<x^2>^-1 * (1-1/n)^2")
print("The resulting contact quotient is a LOWER BOUND, not a measured bond energy.")
print()

central = {}
print("cluster   lower-bound quotient (native)       [radius/gamma band]       record orientation")
for cluster in CLUSTERS:
    label, _, _, _, radii, _, _ = cluster
    value = lower_bound_quotient(cluster, radii[1], GAMMA_CENTRAL)
    band = [lower_bound_quotient(cluster, r, g) for r in radii for g in GAMMA_BAND]
    central[label] = value
    assert abs(value - EXPECTED_ZHZ[label]) < 5e-4, (label, value)
    print(
        f"{label:7s}   {value:6.3f} ZHz"
        f"             [{min(band):6.3f}, {max(band):6.3f}] ZHz"
        f"     (record {record_mev(value):5.2f} MeV)"
    )

values = list(central.values())
spread = max(values) - min(values)
ratio = max(values) / min(values)
print()
print("FROZEN SWING-6 GATES, rerun without moving their thresholds")
gate_a = all(13.0 * MEV_TO_ZHZ <= v <= 19.0 * MEV_TO_ZHZ for v in values)
gate_b = central["alpha"] - central["d"] >= 1.5 * MEV_TO_ZHZ
gate_c = 16.0 * MEV_TO_ZHZ <= central["alpha"] <= 20.0 * MEV_TO_ZHZ
print(f"  T-C6a all central quotients in [3.143,4.594] ZHz (record [13,19] MeV): {'PASS' if gate_a else 'FAIL'}")
print(f"  T-C6b alpha-d enhancement >=0.363 ZHz (record 1.5 MeV):              {'PASS' if gate_b else 'FAIL'}")
print(f"  T-C6c alpha quotient in [3.869,4.836] ZHz (record [16,20] MeV):       {'PASS' if gate_c else 'FAIL'}")
print(f"  central spread: {spread:.3f} ZHz (record {record_mev(spread):.2f} MeV); max/min={ratio:.3f}")
assert not gate_a and gate_b and not gate_c

# The three raw bookkeeping quotients discussed in section 3.
net_pair_quotients = (2.2245662, 8.4817962 / 3.0, 28.2956624 / 6.0)
assert all(a < b for a, b in zip(net_pair_quotients, net_pair_quotients[1:]))
print("  raw B/pair-count sequence is monotone: " +
      " < ".join(f"{x * MEV_TO_ZHZ:.3f} ZHz" for x in net_pair_quotients) +
      " (bookkeeping quotients, not observables)")

# Peak bookkeeping corrections, independent of the cluster-anchor failure.
A_V = 15.75 * MEV_TO_ZHZ
A_SYM_ILLUSTRATIVE = 24.5 * MEV_TO_ZHZ


def valley_z(A: float) -> float:
    p = A_C * A ** (-4.0 / 3.0)
    return (4.0 * A_SYM_ILLUSTRATIVE / A + p) / (
        8.0 * A_SYM_ILLUSTRATIVE / A**2 + 2.0 * p
    )


def valley_b(A: int, ratio: float) -> float:
    z = valley_z(A)
    return (
        A_V
        - ratio * A_V * A ** (-1.0 / 3.0)
        - A_C * z * (z - 1.0) * A ** (-4.0 / 3.0)
        - A_SYM_ILLUSTRATIVE * (A - 2.0 * z) ** 2 / A**2
    )


print()
print("CONDITIONAL PEAK ARITHMETIC (separate from the failed cluster anchor)")
for ratio in (1.13, 1.50):
    balanced = 2.0 * ratio * A_V / A_C - 1.0  # exact Z(Z-1) bookkeeping
    optimized = max(range(4, 271), key=lambda A: valley_b(A, ratio))
    print(
        f"  a_s/a_v={ratio:.2f}: balanced stationary A={balanced:.2f}; "
        f"optimized illustrative valley peak A={optimized}"
    )
assert max(range(4, 271), key=lambda A: valley_b(A, 1.13)) == 60
assert max(range(4, 271), key=lambda A: valley_b(A, 1.50)) == 84

print()
print("VERDICT")
print("  finite-n correction: PASS")
print("  former flat light-cluster anchor: FAIL")
print("  affected gamma/capacity/contact/symmetry chain: QUARANTINE")
print("  raw AME data, unit conversion, pair counting, conditional peak existence: SURVIVE")
