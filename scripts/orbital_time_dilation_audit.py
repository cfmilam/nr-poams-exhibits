#!/usr/bin/env python3
"""Independent math/code regression for orbital-time-dilation.html."""

from __future__ import annotations

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HTML = (ROOT / "orbital-time-dilation.html").read_text(encoding="utf-8")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


# Circular panel: rho-domain must keep both displayed roots real.
for half_step in range(8, 81):
    rho = half_step / 2
    for coefficient in (1, 3):
        radicand = 1 - coefficient / rho
        require(radicand > 0, f"circular radicand failed: rho={rho}, k={coefficient}")
        require(math.isfinite(math.sqrt(radicand)), "circular factor is not finite")

# Given the two printed speed relations, equality is algebraically unique at a=2.
for tenth in range(0, 41):
    a = tenth / 10
    match = math.isclose(a / 2, 1.0, abs_tol=1e-15)
    require(match == math.isclose(a, 2.0, abs_tol=1e-15), "a=2 matching gate failed")

# Ellipse panel: fixed mu, physical a control, real square root throughout the control grid.
MU = 0.002


def solve_kepler(mean_anomaly: float, eccentricity: float) -> float:
    eccentric_anomaly = mean_anomaly
    for _ in range(30):
        eccentric_anomaly -= (
            eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) - mean_anomaly
        ) / (1 - eccentricity * math.cos(eccentric_anomaly))
    return eccentric_anomaly


def true_anomaly(eccentric_anomaly: float, eccentricity: float) -> float:
    return 2 * math.atan2(
        math.sqrt(1 + eccentricity) * math.sin(eccentric_anomaly / 2),
        math.sqrt(1 - eccentricity) * math.cos(eccentric_anomaly / 2),
    )


def ellipse_factor(a: float, eccentricity: float, theta: float) -> tuple[float, float]:
    radius = a * (1 - eccentricity**2) / (1 + eccentricity * math.cos(theta))
    speed_sq = MU * (2 / radius - 1 / a)
    radicand = 1 - 2 * MU / radius + MU / a
    require(math.isclose(radicand, 1 - speed_sq, rel_tol=0, abs_tol=2e-15), "Eq 2.12 substitution failed")
    return math.sqrt(radicand), radicand


minimum_radicand = 1.0
for a_step in range(50, 201, 5):
    a = a_step / 100
    for e_step in range(0, 91):
        eccentricity = e_step / 100
        for theta_step in range(0, 361, 5):
            factor, radicand = ellipse_factor(a, eccentricity, math.radians(theta_step))
            minimum_radicand = min(minimum_radicand, radicand)
            require(radicand > 0 and math.isfinite(factor), "ellipse domain failure")
require(minimum_radicand > 0.92, f"unexpected ellipse margin: {minimum_radicand}")
near = ellipse_factor(0.5, 0.6, 0)[0]
far = ellipse_factor(2.0, 0.6, 0)[0]
require(far > near, "semi-major-axis control does not weaken the displayed deficit")
period_ratio = math.sqrt(2.0**3 / 0.5**3)
require(math.isclose(period_ratio, 8.0, rel_tol=1e-15), "Kepler period scaling failed")


def per_orbit_delta(a: float, eccentricity: float, samples: int = 20000) -> float:
    mean_motion = math.sqrt(MU / a**3)
    d_mean = 2 * math.pi / samples
    coordinate = proper = 0.0
    for i in range(samples):
        mean = (i + 0.5) * d_mean
        ecc_anomaly = solve_kepler(mean, eccentricity)
        theta = true_anomaly(ecc_anomaly, eccentricity)
        factor, _ = ellipse_factor(a, eccentricity, theta)
        dt = d_mean / mean_motion
        coordinate += dt
        proper += factor * dt
    return coordinate - proper


delta_half = per_orbit_delta(0.5, 0.6)
delta_two = per_orbit_delta(2.0, 0.6)
require(delta_half > 0 and delta_two > 0 and not math.isclose(delta_half, delta_two), "ellipse a-control is inert")

# Mercury source-equation check.
G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30
MERCURY_A = 57.90905e9
MERCURY_E = 0.205630
MERCURY_PERIOD_DAYS = 87.9691
ARCSEC_PER_RAD = 206_264.80624709636
mercury_epsilon = G * M_SUN / (MERCURY_A * (1 - MERCURY_E**2) * C**2)
mercury_arcsec_orbit = 6 * math.pi * mercury_epsilon * ARCSEC_PER_RAD
mercury_arcsec_century = mercury_arcsec_orbit * (36_525 / MERCURY_PERIOD_DAYS)
require(abs(mercury_arcsec_orbit - 0.10352) < 5e-5, "Mercury per-orbit value failed")
require(abs(mercury_arcsec_century - 42.982) < 0.01, "Mercury per-century value failed")

# GPS component ledger: altitude plus motion must equal the full relative rate.
M_EARTH = 5.972e24
R_EARTH = 6.371e6
R_GPS = 26.571e6
DAY = 86_400.0
rate_ground = math.sqrt(1 - 2 * G * M_EARTH / (R_EARTH * C**2))
rate_gravity = math.sqrt(1 - 2 * G * M_EARTH / (R_GPS * C**2))
rate_full = math.sqrt(1 - 3 * G * M_EARTH / (R_GPS * C**2))
gravity_us_day = (rate_gravity - rate_ground) * DAY * 1e6
motion_us_day = (rate_full - rate_gravity) * DAY * 1e6
net_us_day = (rate_full - rate_ground) * DAY * 1e6
require(gravity_us_day > 0 and motion_us_day < 0 and net_us_day > 0, "GPS component signs failed")
require(abs(gravity_us_day + motion_us_day - net_us_day) < 1e-9, "GPS ledger does not close")
require(abs(net_us_day - 38.5125) < 0.002, "GPS net failed")

# Source guards for the repaired controls and labels.
required = [
    "Normalized orbital radius",
    "algebraic match between the paper's two stated speed relations",
    "it does not independently derive Eq 3.5 or simulate off-match dynamics",
    "const mu = 0.002",
    "Eq 2.12 with K/c² = -mu/a",
    "const n = Math.sqrt(mu / Math.pow(aModel, 3))",
    "function computeOrbitDelta(a, ecc, samples=1440)",
    "const drift_kinematic = rate_sat_full - rate_sat_grav_only",
    "GPS component ledger failed to close",
    "Source-equation check",
]
for token in required:
    require(token in HTML, f"missing source guard: {token}")
for forbidden in (
    "pretend units",
    "Math.sqrt(val) : 0.001",
    "const drift = (a - 2)",
    "Velocity-only (gravitational) dilation",
    "precession strength <span",
):
    require(forbidden not in HTML, f"obsolete source remains: {forbidden}")

print(f"circular domain: rho=4..40, k=1/3 PASS")
print(f"ellipse minimum radicand: {minimum_radicand:.12f}")
print(f"ellipse per-orbit deficits: a=0.5 {delta_half:.9f}; a=2.0 {delta_two:.9f}")
print(f"Mercury: epsilon={mercury_epsilon:.12e}; {mercury_arcsec_orbit:.6f} arcsec/orbit; {mercury_arcsec_century:.6f} arcsec/century")
print(f"GPS: altitude={gravity_us_day:+.6f}; motion={motion_us_day:+.6f}; net={net_us_day:+.6f} us/day")
print("orbital-time-dilation math/code audit PASS")
