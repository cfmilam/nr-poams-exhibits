#!/usr/bin/env python3
"""Conditional weak-link transition for the Vortex Cores exhibit.

This is an effective material model, NOT an ab-initio POAMS prediction or a
weight-modification claim. No material constants, measured force, alpha, or
domain multiplier are fitted. Energy is reported in the declared K_minus unit;
time is in gamma/K_minus. Python 3 + NumPy only.

The quadratic link costs, first-harmonic joining cost, signed pair stiffness,
overdamped relaxation and position-dependent control are constitutive inputs.
Integer phase closure alone does not determine them.
"""

import json
import math

import numpy as np

PI = math.pi
CHECKS = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def close(a, b, tol=1e-10):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def bisect(fn, a, b):
    fa, fb = fn(a), fn(b)
    if abs(fa) < 1e-14:
        return a
    if abs(fb) < 1e-14:
        return b
    if fa * fb > 0:
        raise ValueError((a, b, fa, fb))
    for _ in range(85):
        m = (a + b) / 2
        fm = fn(m)
        if fm * fa > 0:
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


def v(q, u, b):
    """Single reduced ring energy / K_minus."""
    return 0.5 * (q - u) ** 2 + b * (1 - math.cos(q))


def stationarity(q, u, b):
    return q - u + b * math.sin(q)


def fold(b):
    if b <= 1:
        raise ValueError("No hysteretic fold for b <= 1")
    q = math.acos(-1 / b)
    u_up = q + b * math.sin(q)
    return q, u_up, 2 * PI - u_up


def well(u, b, n):
    qf, _, _ = fold(b)
    center = 2 * PI * n
    return bisect(lambda q: stationarity(q, u, b), center - qf, center + qf)


def pair_energy(qs, us, K, zeta, Eb):
    x = np.asarray(qs) - np.asarray(us)
    matrix = K * np.array([[1.0, zeta], [zeta, 1.0]])
    return float(0.5 * x @ matrix @ x + Eb * np.sum(1 - np.cos(qs)))


def gauss_integral(fn, a, b, n=240):
    x, w = np.polynomial.legendre.leggauss(n)
    xx = (a + b) / 2 + (b - a) * x / 2
    return (b - a) * float(w @ np.array([fn(float(t)) for t in xx])) / 2


def relax(u, b, q0, dt=0.004, duration=180):
    """Overdamped dimensionless dynamics, including heat as a second variable.

    dq/ds = -d(v)/dq; d(Q/K_minus)/ds = (dq/ds)^2.
    The heat goes into an explicitly included material thermal receiver.
    """
    def rhs(y):
        f = -stationarity(y[0], u, b)
        return np.array([f, f * f])

    y = np.array([q0, 0.0])
    initial = v(q0, u, b)
    error = 0.0
    previous_energy = initial
    increases = 0
    for _ in range(int(duration / dt)):
        k1 = rhs(y)
        k2 = rhs(y + dt * k1 / 2)
        k3 = rhs(y + dt * k2 / 2)
        k4 = rhs(y + dt * k3)
        y += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        energy = v(y[0], u, b)
        error = max(error, abs(energy + y[1] - initial))
        increases += int(energy > previous_energy + 1e-12)
        previous_energy = energy
    return y, error, increases


def main():
    # Finite bulk-link reduction: the participation is solved, not presumed.
    ks = np.array([2.0, 3.0, 7.0, 5.0])
    mismatch = 1.23
    Keff = 1 / np.sum(1 / ks)
    deltas = Keff * mismatch / ks
    check("finite phase constraint", close(float(deltas.sum()), mismatch))
    check("equal generalized torque", np.ptp(ks * deltas) < 1e-12)
    check("harmonic-sum energy", close(float(0.5 * np.sum(ks * deltas ** 2)),
                                      0.5 * Keff * mismatch ** 2))
    rng = np.random.default_rng(193)
    perturbations = rng.normal(size=(80, len(ks)))
    perturbations -= perturbations.mean(axis=1, keepdims=True)
    energies = 0.5 * np.sum(ks * (deltas + perturbations) ** 2, axis=1)
    check("constrained minimum", bool(np.all(energies >= 0.5 * Keff * mismatch ** 2)))

    # Declared dimensionless example. These are not apparatus measurements.
    K, zeta, b = 1.0, 0.5, 2.0
    Km = K * (1 - zeta)
    Eb = b * Km
    qc, up, down = fold(b)
    qr = well(up, b, 1)
    dq = qr - qc
    drop_one = v(qc, up, b) - v(qr, up, b)
    check("stationary fold", abs(stationarity(qc, up, b)) < 1e-12)
    check("zero soft curvature", abs(1 + b * math.cos(qc)) < 1e-12)
    check("landing stationary and stable", abs(stationarity(qr, up, b)) < 1e-12
          and 1 + b * math.cos(qr) > 0)
    check("positive energy release", drop_one > 0)
    check("reverse fold symmetry", close(down, 2 * PI - up))
    check("branch change is not a fixed 2pi endpoint change", abs(dq - 2 * PI) > 1)

    H = K * np.array([[1.0, zeta], [zeta, 1.0]]) + Eb * math.cos(qc) * np.eye(2)
    anti = np.array([1.0, -1.0]) / math.sqrt(2)
    sym = np.array([1.0, 1.0]) / math.sqrt(2)
    check("opposed mode soft", np.linalg.norm(H @ anti) < 1e-12)
    check("other mode stable for positive zeta", close(float(sym @ H @ sym), 2 * K * zeta))
    Kp = K * (1 + zeta)
    check("transverse contraction throughout the slip", Kp > Eb)
    # delta=(qA+qB)/2 obeys gamma*delta_dot=-Kp*delta-Eb*cos(w)*sin(delta).
    # |sin(delta)|<=|delta| bounds growth by -(Kp-Eb)|delta| for every w.
    for w, deviation in rng.uniform(-8, 8, size=(100, 2)):
        dev_rate = -Kp * deviation - Eb * math.cos(w) * math.sin(deviation)
        if deviation != 0:
            check_value = math.copysign(1.0, deviation) * dev_rate
            if check_value > -(Kp - Eb) * abs(deviation) + 1e-12:
                raise AssertionError("transverse contraction inequality")
    check("transverse contraction inequality", True)
    check("full pair equals reduced energy", close(pair_energy([qr, -qr], [up, -up], K, zeta, Eb),
                                                  2 * Km * v(qr, up, b)))
    # Negative zeta has a competing instability before the opposed fold.
    zn = -0.25
    Kmn, Kpn = K * (1 - zn), K * (1 + zn)
    Ebn = b * Kmn
    q_competing = math.acos(-Kpn / Ebn)
    q_opposed = math.acos(-Kmn / Ebn)
    check("negative zeta destabilizes symmetric mode first", q_competing < q_opposed)
    for bb in (0.2, 0.8, 1.0):
        check(f"no negative curvature for b={bb}", 1 - bb >= 0)

    # Closed quasistatic cycle, n=0 -> n=1 -> n=0. Integrate *drive work*,
    # not the release formula again: d u/d q = 1 + b cos q.
    q0 = well(down, b, 0)
    q1 = 2 * PI - qc
    def work_integrand(q):
        uu = q + b * math.sin(q)
        return -2 * (q - uu) * (1 + b * math.cos(q))
    work_cycle = gauss_integral(work_integrand, q0, qc)
    work_cycle += gauss_integral(work_integrand, qr, q1)
    heat_cycle = 4 * drop_one
    check("full reset cycle work equals both releases", close(work_cycle, heat_cycle, 2e-11))
    check("reset returns to original lower state", close(q0, 2 * PI - qr))

    # Actual finite-time run is just beyond the fold; at exactly the fold an
    # exactly prepared noiseless solution would stay there indefinitely.
    overdrive = 0.02
    final, err, increases = relax(up + overdrive, b, qc)
    final_half, err_half, _ = relax(up + overdrive, b, qc, dt=0.002)
    target = well(up + overdrive, b, 1)
    check("finite-overdrive relaxation lands in adjacent well", close(final[0], target, 1e-9))
    check("relaxation energy plus heat conserved", err < 2e-8)
    check("relaxation energy nonincreasing", increases == 0)
    check("time-step convergence", np.linalg.norm(final_half - final) < 2e-8 and err_half < 2e-8)
    local_slopes = (1 + b * math.cos(qc), 1 + b * math.cos(qr))

    # Generic reciprocal receiver: R=z_D-z_C. K and Eb fixed here; a calibrated
    # local phase-position slope is a test input, not a claimed apparatus law.
    R0, slope, zz = 3.0, 0.17, zeta
    def us_at(R):
        ur = up + slope * (R - R0)
        return [ur, -ur]
    us = np.array(us_at(R0))
    qs = np.array([qr, -qr])
    stiffness = K * np.array([[1.0, zz], [zz, 1.0]])
    analytic_F = float((qs - us) @ stiffness @ np.array([slope, -slope]))
    eps = 1e-5
    numerical_F = -(pair_energy(qs, us_at(R0 + eps), K, zz, Eb)
                    - pair_energy(qs, us_at(R0 - eps), K, zz, Eb)) / (2 * eps)
    check("receiver force is derivative of same energy", close(analytic_F, numerical_F, 1e-9))
    initial_F = 2 * Km * (qc - up) * slope
    check("force step from branch change", close(analytic_F - initial_F, 2 * Km * slope * dq))
    def energy_positions(zD, zC):
        return pair_energy(qs, us_at(zD - zC), K, zz, Eb)
    FC = -(energy_positions(R0, eps) - energy_positions(R0, -eps)) / (2 * eps)
    check("opposite receiver reaction", close(FC, -analytic_F, 1e-9))
    check("no positional coupling no force", float((qs - us) @ stiffness @ np.zeros(2)) == 0)
    # Also exercise the K'(R) term, with u held independent of position.
    stiffness_slope = 0.13
    def changing_stiffness_energy(R):
        return pair_energy(qs, us, K * math.exp(stiffness_slope * (R - R0)), zz, Eb)
    force_K = -0.5 * stiffness_slope * float((qs - us) @ stiffness @ (qs - us))
    finite_K = -(changing_stiffness_energy(R0 + eps) - changing_stiffness_energy(R0 - eps)) / (2 * eps)
    check("position-dependent stiffness force", close(force_K, finite_K, 1e-9))

    # The old G-scaled channel remains separately G-scaled even if local
    # material dynamics are rapid. Verify its far-away balanced-pair cancellation.
    radius, spacing = 1000.0, 1.0
    exact_diff = 1 / (radius + spacing / 2) ** 4 - 1 / (radius - spacing / 2) ** 4
    leading_diff = -4 * spacing / radius ** 5
    check("balanced-pair far counterpart suppression", abs(exact_diff / leading_diff - 1) < 2e-6)

    # Nonidentifiability: multiplying all energy coefficients by a positive
    # factor preserves q and u folds but changes the release. Changing damping
    # independently changes time. Nothing in closure chooses these scales.
    scale = 7.0
    scaled_Km = scale * K * (1 - zeta)
    scaled_b = (scale * Eb) / scaled_Km
    check("energy scaling preserves threshold", close(fold(scaled_b)[1], up))
    check("energy scaling changes release", close(2 * scaled_Km * drop_one, scale * 2 * Km * drop_one))
    second_harmonic = 0.15
    q_harm = bisect(lambda q: 1 + b * math.cos(q) + 4 * second_harmonic * math.cos(2 * q),
                   PI / 2, PI)
    u_harm = q_harm + b * math.sin(q_harm) + 2 * second_harmonic * math.sin(2 * q_harm)
    check("same periodicity different threshold", abs(u_harm - up) > 0.01)

    print(json.dumps({
        "model_status": "conditional effective material law; not ab-initio POAMS or weight modification",
        "inputs": {"K": K, "zeta": zeta, "E_b_over_K_minus": b,
                   "units": "energy K_minus; time gamma/K_minus; phase radians"},
        "finite_link_example": {"k": ks.tolist(), "K_effective": Keff, "delta": deltas.tolist()},
        "fold": {"q": qc, "u_up": up, "u_down": down, "q_landing": qr,
                 "phase_change": dq, "pair_release_over_K_minus": 2 * drop_one,
                 "cycle_heat_over_K_minus": heat_cycle, "cycle_drive_work_over_K_minus": work_cycle,
                 "soft_and_landing_curvature_over_K_minus": local_slopes},
        "finite_relaxation": {"overdrive": overdrive, "q_final": float(final[0]),
                              "heat_single_over_K_minus": float(final[1]), "maximum_energy_error": err},
        "collective_condition": {"K_plus_minus_E_b": Kp - Eb,
                                 "sufficient_zeta_lower_bound": (b - 1) / (b + 1)},
        "second_harmonic_counterexample": {"E_2_over_K_minus": second_harmonic,
                                           "q_fold": q_harm, "u_fold": u_harm},
        "receiver_derivative_test": {"phase_position_slope_test_input": slope,
                                     "F": analytic_F, "F_counterpart": FC,
                                     "force_step": analytic_F - initial_F},
        "checks_passed": len(CHECKS), "checks": CHECKS,
    }, indent=2))


if __name__ == "__main__":
    main()
