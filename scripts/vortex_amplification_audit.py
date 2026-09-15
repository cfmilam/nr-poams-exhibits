#!/usr/bin/env python3
"""Vortex Cores: finite-ledger checks, not a fitted apparatus simulation.

Run: python3 scripts/vortex_amplification_audit.py
All numerical cases below are declared examples. No alpha value, measured
anomalous force, domain size, or transition rate is inferred from a target.
"""

import json
import math
from fractions import Fraction as Q

H = 6.62607015e-34
C = 299792458.0
HBAR = H / (2 * math.pi)
G = 6.67430e-11
G0 = 9.80665
EARTH_MASS = 5.9722e24
EARTH_RADIUS = 6.371e6


def close(a, b, tolerance=1e-12):
    assert math.isclose(a, b, rel_tol=tolerance, abs_tol=0), (a, b)


def branch_energy(k, u, spacing, scale):
    """Counterexample family, NOT a proposed device law; all inputs rational."""
    return scale * (k - u / spacing) ** 2


def participation(weights):
    """Descriptive participation of nonnegative event weights, not phase order."""
    return sum(weights) ** 2 / sum(w * w for w in weights)


def recoil(event_joules, rate_hz, eta):
    """ONLY the corpus's one-way p=E/c exchange channel, with signed eta."""
    assert event_joules >= 0 and rate_hz >= 0 and abs(eta) <= 1
    power = event_joules * rate_hz
    return {"power_W": power, "force_N": eta * power / C}


def main():
    report = {"scope": "Conditional ledger identities and counterexamples; no device prediction."}

    # Same integer states, arbitrarily different crossing locations and gaps.
    # Positive rate E/hbar represents each positive energy in the stated rate ledger.
    threshold_examples = []
    for spacing in [Q(1), Q(7), Q(100)]:
        for scale in [Q(1), Q(5)]:
            u = spacing / 2
            assert branch_energy(0, u, spacing, scale) == branch_energy(1, u, spacing, scale)
            assert branch_energy(0, u / 2, spacing, scale) < branch_energy(1, u / 2, spacing, scale)
            assert branch_energy(0, u * Q(3, 2), spacing, scale) > branch_energy(1, u * Q(3, 2), spacing, scale)
            threshold_examples.append({"spacing": float(spacing), "scale": float(scale), "crossing": float(u)})
    report["threshold_nonuniqueness"] = threshold_examples

    # Finite difference of sum_i n_i * omega_i; fixed frequency is a special case.
    n, dn = [2, 3], [1, -1]
    w, dw = [Q(5), Q(7)], [Q(1, 3), Q(-1, 2)]
    direct = sum((ni + di) * (wi + dwi) - ni * wi for ni, di, wi, dwi in zip(n, dn, w, dw))
    split = sum(wi * di + (ni + di) * dwi for ni, di, wi, dwi in zip(n, dn, w, dw))
    assert direct == split
    report["finite_rate_ledger_identity"] = str(direct)

    # Exactly the same phase order (all phases zero) permits distinct event supports.
    size = 64
    phase_order = abs(sum(complex(1, 0) for _ in range(size))) ** 2 / size ** 2
    assert phase_order == 1
    one, all_turns = [Q(1)] + [Q(0)] * (size - 1), [Q(1)] * size
    assert participation(one) == 1 and participation(all_turns) == size
    # A synchronized exchange between equal-rate entries can cancel its net energy.
    assert sum([1, -1]) == 0
    report["coherence_does_not_fix_participation"] = {"phase_order": phase_order, "one_change": int(participation(one)), "all_change": int(participation(all_turns))}

    # A coarse-grained vector dot product already contains each microscopic term.
    vectors = [(Q(1), Q(0), Q(0)), (Q(-1), Q(2), Q(0)), (Q(0), Q(0), Q(3))]
    partner = (Q(2), Q(3), Q(5))
    total = tuple(sum(v[i] for v in vectors) for i in range(3))
    macro = sum(a * b for a, b in zip(total, partner))
    micro = sum(sum(a * b for a, b in zip(v, partner)) for v in vectors)
    assert macro == micro
    report["no_second_domain_multiplier"] = str(macro)

    # The carried -2 dot-product/r^3 kernel is outward for opposed axes.
    r, coefficient = Q(2), Q(3)
    radial_aligned = -3 * coefficient / r ** 4
    radial_opposed = 3 * coefficient / r ** 4
    assert radial_aligned < 0 < radial_opposed
    report["pair_kernel_radial_sign"] = {"aligned": "inward", "opposed": "outward"}

    # Rotational scalar dot product alone is not the only invariant:
    # at L1 perpendicular L2 their projections on rhat can still be nonzero.
    L1_dot_L2 = 0.0
    projection_product = (1 / math.sqrt(2)) ** 2
    assert L1_dot_L2 == 0 and projection_product > 0
    report["independent_rotational_invariant_example"] = {"dot_product": L1_dot_L2, "radial_projection_product": projection_product}

    # Internal pair impulses cancel; an external counterpart is explicit.
    endpoint_impulses = [Q(2), Q(-2), Q(7, 3), Q(-7, 3)]
    assert sum(endpoint_impulses) == 0
    # Reflection-symmetric directional exchange: equal upward/downward energy.
    assert (Q(5) - Q(5)) == 0
    report["balanced_exchange_net_impulse"] = 0

    # Fixed output power: event size changes rate, not mean recoil.
    fixed_power, eta, nu = 6.0, 0.5, 1e10
    events = []
    for q in [1, 1e6, 1e14, 1e20]:
        energy = q * H * nu
        rate = fixed_power / energy
        out = recoil(energy, rate, eta)
        close(out["force_N"], eta * fixed_power / C)
        events.append({"equal_gap_steps": q, "event_J": energy, "event_rate_Hz": rate, **out})
    report["fixed_power_cases"] = events

    # Reproduce the old rounded event assumption before any 2*pi convention fix.
    old_rows = []
    for rate, printed in [(1e3, 1.7e-13), (1e6, 1.7e-10), (1e9, 1.7e-7)]:
        out = recoil(1e-4, rate, 0.5)
        assert 970 < out["force_N"] / printed < 1000
        old_rows.append({"rate_Hz": rate, "old_printed_force_N": printed, **out})
    report["old_rounded_table_recalculated"] = old_rows
    angular_event = 1e20 * HBAR * 1e10
    cycle_event = 1e20 * H * 1e10
    close(cycle_event / angular_event, 2 * math.pi)
    report["frequency_conventions"] = {"omega_1e10_rad_s_event_J": angular_event, "nu_10_GHz_event_J": cycle_event}
    report["ten_GHz_conditional_table"] = [{"rate_Hz": rate, **recoil(cycle_event, rate, 0.5)} for rate in [1e3, 1e6, 1e9]]

    # References are illustrative readout targets, NOT universal torsion sensitivity.
    tiny_force = 1e-9 * 1e-3 * G0  # 10^-9 gram mass equivalent -> N
    report["reference_targets"] = [{"force_N": force, "required_output_power_W_at_eta_half": force * C / 0.5} for force in [1e-8, tiny_force]]
    report["old_threshold_ratio"] = 1e-8 / tiny_force

    # Preserve prior continuous-model arithmetic, explicitly holding L vectors fixed
    # in the force derivative and then evaluating on a circular reference orbit.
    mass = 1.0
    L_orbit = mass * math.sqrt(G * EARTH_MASS * EARTH_RADIUS)
    factor = 6 * G * mass / (C ** 2 * EARTH_RADIUS)
    continuous = []
    for L in [1e3, 1e5, 1e7]:
        beta = L / L_orbit
        fraction = beta * factor
        direct = 6 * L * L_orbit / (C ** 2 * EARTH_RADIUS ** 2 * mass * EARTH_MASS)
        close(fraction, direct)
        continuous.append({"assumed_L_J_s": L, "beta": beta, "abs_delta_g_over_g": fraction, "mass_equivalent_g": fraction * mass * 1000})
    report["continuous_model_fixed_L_derivative"] = continuous
    # Substituting L_orbit(r) before differentiating changes coefficient 6 to 5;
    # neither calculation derives a physically selected constrained trajectory.
    report["circular_substitution_derivative_ratio"] = 5 / 6

    # Burst detection needs mechanical bandwidth. Example oscillator response ratio.
    # H(f)/H(0) = 1 / [1-(f/f0)^2 + i*f/(Q*f0)] (calibrated readout example only).
    f0, quality, f_event = 1.0, 100.0, 1e9
    response = 1 / abs(complex(1 - (f_event / f0) ** 2, f_event / (quality * f0)))
    close(response, 1e-18)
    report["illustrative_1Hz_readout_at_1GHz_relative_response"] = response
    report["result"] = "PASS: all stated identities, arithmetic and counterexamples checked"
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
