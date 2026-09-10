#!/usr/bin/env python3
"""Deterministic mathematics-and-logic checks for The Patient Sun."""

from __future__ import annotations

import json
import math
from pathlib import Path

L_SUN = 3.828e26
M_SUN = 1.98847e30
R_SUN = 6.957e8
C = 299_792_458.0
YEAR = 365.25 * 86400
H = 6.62607015e-34
K_B = 1.380649e-23


def main() -> None:
    source = Path(__file__).resolve().parents[1] / "the-patient-sun.html"
    html = source.read_text(encoding="utf-8")

    specific_power = L_SUN / M_SUN
    fractional_rate = L_SUN / (M_SUN * C**2)
    tenure_seconds = 0.007 * 0.1 * M_SUN * C**2 / L_SUN
    processed_mass_rate = L_SUN / (0.007 * C**2)
    published_mass_rate = L_SUN / C**2
    thermal_turn_rate = K_B * 1.57e7 / H
    handoffs = (R_SUN / 1e-3) ** 2
    detention_seconds = R_SUN**2 / (1e-3 * C)
    ballistic_seconds = R_SUN / C
    retail_count = 6.463e21 / 314e12

    assert math.isclose(specific_power, 1.925e-4, rel_tol=5e-4)
    assert math.isclose(fractional_rate, 2.142e-21, rel_tol=5e-4)
    assert math.isclose(tenure_seconds / YEAR / 1e9, 10.37, rel_tol=5e-3)
    assert math.isclose(processed_mass_rate, 6.085e11, rel_tol=5e-3)
    assert math.isclose(published_mass_rate, 4.259e9, rel_tol=5e-3)
    assert math.isclose(thermal_turn_rate / 1e15, 327.17, rel_tol=5e-3)
    assert 4.8e23 < handoffs < 4.9e23
    assert 5.0e4 < detention_seconds / YEAR < 5.2e4
    assert math.isclose(ballistic_seconds, 2.3203, rel_tol=5e-4)
    assert 2.05e7 < retail_count < 2.07e7

    # The displayed toy is self-consistent over its deliberately local slider domain.
    for log_s in (-0.5, 0.0, 0.5):
        s = 10**log_s
        temperature_mk = 15.7 * s ** -0.25
        assert 11.7 < temperature_mk < 21.0
    for log_k in (-0.5, 0.0, 0.5):
        kappa = 10**log_k
        luminosity = L_SUN / kappa
        tenure_gyr = 10.5 * kappa
        assert math.isclose(luminosity * tenure_gyr, L_SUN * 10.5, rel_tol=1e-15)

    required = (
        "dominance is a POAMS hypothesis",
        "Historical motivation, not proof",
        "illustrative local algebra",
        "energy-equivalent average count",
        "within stated experimental and",
        "\\Delta m^2 L/E",
        "neutral-atom mass bookkeeping",
    )
    forbidden = (
        "historical proof that the exit sets the rate",
        "exactly the rate the surface publishes",
        "century-delayed retail",
        "longer than the present age of the universe's stars",
        "Make the entry a million times easier",
        "because it beats in transit",
    )
    assert all(text in html for text in required)
    assert all(text not in html for text in forbidden)

    print(json.dumps({
        "status": "PASS",
        "specific_power_W_per_kg": specific_power,
        "fractional_rate_per_s": fractional_rate,
        "fuel_budget_Gyr": tenure_seconds / YEAR / 1e9,
        "processed_mass_kg_per_s": processed_mass_rate,
        "published_mass_kg_per_s": published_mass_rate,
        "thermal_turn_rate_PHz": thermal_turn_rate / 1e15,
        "constant_mfp_handoffs": handoffs,
        "constant_mfp_detention_years": detention_seconds / YEAR,
        "ballistic_seconds": ballistic_seconds,
        "retail_energy_equivalent_count": retail_count,
    }, indent=2))


if __name__ == "__main__":
    main()
