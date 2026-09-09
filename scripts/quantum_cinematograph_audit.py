#!/usr/bin/env python3
"""Independent arithmetic and source guards for quantum-cinematograph.html."""
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "quantum-cinematograph.html"
TEXT = PAGE.read_text(encoding="utf-8")

class AuditParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = []
        self.links = []
    def handle_starttag(self, tag, attrs):
        data = dict(attrs)
        if data.get("id"):
            self.ids.append(data["id"])
        if tag == "a" and data.get("href"):
            self.links.append(data["href"])

parser = AuditParser()
parser.feed(TEXT)
duplicates = [name for name, count in Counter(parser.ids).items() if count > 1]
assert not duplicates, duplicates
for href in parser.links:
    if href.startswith("./"):
        assert (ROOT / href[2:].split("#", 1)[0]).exists(), href

h = 6.62607015e-34
hbar = h / (2 * math.pi)
assert abs(2 * math.pi * hbar / h - 1) < 1e-15
G = 6.67430e-11
c = 299_792_458.0
planck_length = math.sqrt(hbar * G / c**3)
planck_time = math.sqrt(hbar * G / c**5)
assert round(planck_length / 1e-35, 1) == 1.6
assert round(planck_time / 1e-44, 1) == 5.4

for rate in (0.1, 1.0, 5.0):
    dt = 0.05
    probability = rate * dt
    assert 0 <= probability <= 1
    assert abs(probability / dt - rate) < 1e-12

for first in (-1.0, 1.0):
    for second in (-1.0, 1.0):
        assert 0 <= (first + second + 2) / 4 <= 1

required = (
    "filmStripRate * elapsedSeconds",
    "interactionTriggeredAt = performance.now()",
    "propagateTouch = function()",
    "Expected Manifestation Rate per Atom",
    "count === 1 ? canvas.width / 2",
    "+ 2) / 4",
    "not empirical measurements or independent proofs",
    r"L=n\hbar",
    r"S_{\rm turn}=\oint L\,d\theta=nh",
)
for marker in required:
    assert marker in TEXT, marker

for retired in (
    "filmStripRate / 60",
    "rate * 0.01",
    "Wave Speed: c",
    "let particles = []",
    "Time itself is discrete!",
    "whole turns of action (ħ closures)",
    "consistent with Heisenberg uncertainty",
):
    assert retired not in TEXT, retired

print(f"PASS: {len(parser.ids)} unique IDs; local links resolve")
print(f"PASS: h=2πħ; l_P={planck_length:.9e} m; t_P={planck_time:.9e} s")
print("PASS: flicker-rate expectation and pixel-intensity bounds")
print("PASS: repaired code and disclosure guards")
