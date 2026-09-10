#!/usr/bin/env python3
"""Deterministic mathematics-and-logic release checks for poams-electricity.html."""
from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "poams-electricity.html"
text = PAGE.read_text(encoding="utf-8")

required = [
    "The record's claim (i)",
    "replication required",
    "In a finite instrument, “none” means $|I_A|&lt;\\delta I$",
    "Conceptual ledger animation, not an experimental trace or circuit solver",
    "P=\\sum_k\\dot N_k h\\nu_k",
    "idealized Planck-spectrum equivalent",
    "one completed return relation",
    "That numerical compatibility is necessary, but it is not independent evidence",
    "pendingDelay=lenMeters()/299792458",
    "requestAnimationFrame(tick)",
]
for token in required:
    assert token in text, f"missing required disclosure or guard: {token}"

forbidden = [
    "about ten atomic spacings before handing off",
    "outside the conductor entirely",
    "Nothing traveled, nothing radiated, nothing moved but a needle",
    '"satellite hop"',
    "300,000 km</span>",
    "the same iron-garnet class POAMS selects from first principles",
    "No such transient exists.",
    "not one measured number changed",
    "tick();",
]
for token in forbidden:
    assert token not in text, f"retired overclaim or defect remains: {token}"

# Independent arithmetic on the page's stated inputs.
I = 1.0
n = 8.49e28
A = 1e-6
e = 1.602176634e-19
me = 9.1093837139e-31
L = 1000.0
c = 299_792_458.0
h = 6.62607015e-34
k = 1.380649e-23
T = 2700.0
P = 60.0

vd = I / (n * A * e)
assert math.isclose(vd * 1e3, 0.07352, rel_tol=8e-5)
assert math.isclose(L / vd / 86400, 157.44, rel_tol=8e-5)
drift_energy = 0.5 * me * vd**2 * n * A * L
assert math.isclose(drift_energy, 2.09e-13, rel_tol=0.006)
assert math.isclose(drift_energy / 100 * 1e15, 2.09, rel_tol=0.006)

assert math.isclose(e / h / 1e12, 241.798924, rel_tol=1e-8)
mean_nu = 2.701 * k * T / h
equiv_rate = P / (2.701 * k * T)
assert math.isclose(mean_nu / 1e12, 151.96, rel_tol=8e-5)
assert math.isclose(equiv_rate, 5.96e20, rel_tol=8e-4)
assert math.isclose(equiv_rate / 120, 4.97e18, rel_tol=1e-3)

# The table is explicitly a vacuum baseline; check each printed value.
assert math.isclose(1 / c * 1e9, 3.336, rel_tol=2e-4)
assert math.isclose(300_000 / c * 1e3, 1.001, rel_tol=5e-4)
assert math.isclose(35_786_000 / c * 1e3, 119.37, rel_tol=5e-5)
assert math.isclose(71_572_000 / c * 1e3, 238.74, rel_tol=5e-5)
assert math.isclose(384_400_000 / c, 1.282, rel_tol=2e-4)
assert math.isclose(1.468 * 1000 / c * 1e6, 4.897, rel_tol=1e-4)

# Dimensional closure: count-rate * action * turn-rate = power.
# Base exponents are (kg, m, s).
count_rate = (0, 0, -1)
action = (1, 2, -1)
turn_rate = (0, 0, -1)
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
assert add(count_rate, add(action, turn_rate)) == (1, 2, -3)


class AuditParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: list[str] = []
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        data = dict(attrs)
        if data.get("id"):
            self.ids.append(data["id"])
        if tag == "a" and data.get("href"):
            self.hrefs.append(data["href"])


parser = AuditParser()
parser.feed(text)
assert len(parser.ids) == len(set(parser.ids)), "duplicate HTML id"
for href in parser.hrefs:
    if href.startswith(("./", "../")):
        target = (PAGE.parent / href.split("#", 1)[0]).resolve()
        assert target.exists(), f"missing local link target: {href}"

for op, cl in [("(", ")"), ("{", "}"), ("[", "]")]:
    assert text.count(op) == text.count(cl), f"unbalanced {op}{cl}"

print(
    "PASS electricity: arithmetic, dimensions, evidence labels, latency table, "
    "IDs, local links, and delayed-reading guards"
)
