#!/usr/bin/env python3
"""Deterministic checks for the Heat & Coherence normalized ledger."""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

N = 64
KMAX = 5.5
OMEGA = 1.15
DNOISE = 3.1
POWER = 0.95
DT = 0.008


def order(phases: list[float]) -> tuple[float, float]:
    c = sum(math.cos(x) for x in phases) / N
    s = sum(math.sin(x) for x in phases) / N
    return math.hypot(c, s), math.atan2(s, c)


def run(coherence: float) -> dict[str, float]:
    seed = 42000 + round(coherence * 1000)
    rng = random.Random(seed)
    phases = [rng.random() * 2 * math.pi for _ in range(N)]
    previous_resultant = 0.0
    old_samples: list[float] = []
    new_samples: list[float] = []
    r_samples: list[float] = []

    for step in range(30000):
        r0, psi0 = order(phases)
        ordered_share = coherence * POWER
        churn_share = (1 - coherence) * POWER
        directed = OMEGA * ordered_share
        noise_amplitude = math.sqrt(2 * DNOISE * churn_share * DT)
        increments = []
        for index, phase in enumerate(phases):
            increment = (
                directed * DT
                + KMAX * coherence * r0 * math.sin(psi0 - phase) * DT
                + noise_amplitude * rng.gauss(0, 1)
            )
            increments.append(increment)
            phases[index] += increment

        r1, psi1 = order(phases)
        resultant_increment = psi1 - previous_resultant
        while resultant_increment > math.pi:
            resultant_increment -= 2 * math.pi
        while resultant_increment < -math.pi:
            resultant_increment += 2 * math.pi
        previous_resultant = psi1
        mean_increment = sum(increments) / N

        if step >= 5000:
            old_samples.append(
                sum((x - resultant_increment) ** 2 for x in increments) / N / DT
            )
            new_samples.append(
                sum((x - mean_increment) ** 2 for x in increments) / N / DT
            )
            r_samples.append(r1)

    expected_noise = 2 * (1 - coherence) * DNOISE * POWER
    return {
        "coherence": coherence,
        "mean_r": sum(r_samples) / len(r_samples),
        "old_churn": sum(old_samples) / len(old_samples),
        "repaired_churn": sum(new_samples) / len(new_samples),
        "noise_expectation": expected_noise,
    }


def main() -> None:
    source = Path(__file__).resolve().parents[1] / "heat-and-coherence.html"
    html = source.read_text(encoding="utf-8")

    required = (
        "P<sub>churn</sub> = (1 − c)P",
        "P<sub>ordered</sub> = cP",
        "held fixed; not computed here",
        "not a scale output",
        "meanDphi",
        "not a material relaxation measurement",
    )
    forbidden = ("killTimer", "dphi[m] - d", "weight displaced")
    assert all(text in html for text in required)
    assert all(text not in html for text in forbidden)

    for coherence in (0.0, 0.04, 0.5, 0.98, 1.0):
        for power in (0.3, 0.95, 1.5):
            assert math.isclose(
                (1 - coherence) * power + coherence * power,
                power,
                rel_tol=0,
                abs_tol=1e-15,
            )

    results = [run(c) for c in (0.06, 0.5, 0.75, 0.98)]
    for result in results:
        assert 0 <= result["mean_r"] <= 1
        expected = result["noise_expectation"]
        if expected:
            relative_error = abs(result["repaired_churn"] - expected) / expected
            assert relative_error < 0.06, result
    assert results[0]["old_churn"] > 3 * results[0]["noise_expectation"]
    assert results[0]["repaired_churn"] < 0.35 * results[0]["old_churn"]

    print(json.dumps({"status": "PASS", "runs": results}, indent=2))


if __name__ == "__main__":
    main()
