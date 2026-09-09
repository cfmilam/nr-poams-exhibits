#!/usr/bin/env python3
"""Deterministic release checks for antenna-emission.html."""
from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "antenna-emission.html"
text = PAGE.read_text(encoding="utf-8")

required = [
    "h=2\\pi\\hbar",
    "P=\\Gamma\\Delta E=\\Gamma h\\nu",
    "constitutive calculation of $\\Gamma$ open",
    "Conceptual schematic, not a solver or evidentiary trace",
    "contested research lead",
    "No equivalence to a particular historical direct-action force law is asserted",
    "twist-sense supplies only the handedness component",
    "display pixels/s; no physical scale",
]
for token in required:
    assert token in text, f"missing required disclosure: {token}"

forbidden = [
    "the experiments side with Ampère",
    "The force Maxwell's field cannot hold",
    "Power out = AM leaving the ledger",
    "An antenna receives with exactly the pattern it transmits",
    "later preferred in the direct-action work of Feynman and Wheeler",
    "longitudinal instantaneous hand-off → absorber",
]
for token in forbidden:
    assert token not in text, f"retired overclaim remains: {token}"

# Dimensional ledger in base exponents (kg, m, s).
J = (1, 2, -2)
h = (1, 2, -1)       # J s
nu = (0, 0, -1)      # s^-1
gamma = (0, 0, -1)   # closures s^-1

add = lambda a, b: tuple(x + y for x, y in zip(a, b))
assert add(h, nu) == J
assert add(gamma, add(h, nu)) == (1, 2, -3)  # watt
assert add(h, gamma) == J                    # dL/dt is energy, not power

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

# Every inline script must at least have balanced delimiters; Node performs syntax parsing.
for op, cl in [("(", ")"), ("{", "}"), ("[", "]")]:
    assert text.count(op) == text.count(cl), f"unbalanced {op}{cl}"

# The animation must advance by elapsed seconds, not frames.
assert "function frame(now)" in text
assert "t+=dt" in text
assert "t+=1" not in text

print("PASS antenna-emission: dimensions, disclosures, IDs, local links, and timing guards")
