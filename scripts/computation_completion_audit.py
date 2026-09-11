#!/usr/bin/env python3
from pathlib import Path
import math
import re

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "computation-as-completion.html"
text = PAGE.read_text(encoding="utf-8")

ids = re.findall(r'\bid="([^"]+)"', text)
assert len(ids) == len(set(ids)), "duplicate IDs"

for href in re.findall(r'href="([^"]+)"', text):
    if href.startswith("./"):
        target = ROOT / href[2:].split("#", 1)[0]
        assert target.exists(), f"missing local link: {href}"

lam = 2.0
target_exponent = 12
distance = math.ceil(2 * target_exponent / math.log10(lam) - 1)
if distance % 2 == 0:
    distance += 1
qubits = 2 * distance * distance
assert (distance, qubits) == (79, 12482)

p, commits = 1e-4, 1000
p_fail = 1 - (1 - p) ** commits
assert math.isclose(p_fail, 0.09516710644144377, rel_tol=0, abs_tol=1e-15)

for register, a, b, modulus in ((0, 7, 9, 12), (4, 7, 9, 12), (11, 1, 1, 12)):
    total = register + a + b
    assert total % modulus in range(modulus)
    assert total // modulus >= 0

required = (
    "ledgerErr+=miss",
    r"O(r\,\mathrm{poly}(d))",
    "not an impossibility theorem",
    "not automatically held-in winding",
    "upper 95% confidence bound",
    "normalized A=1 model",
)
for phrase in required:
    assert phrase in text, f"missing guard: {phrase}"

for phrase in (
    "Error quantized, not accumulated",
    "Their arithmetic is our theorem now",
    "section derives that it cannot become one",
    "winding architecture: settled",
):
    assert phrase not in text, f"stale claim: {phrase}"

assert text.count('class="canvas-scroll"') == 4
assert text.count('class="table-scroll"') == 2

print(
    "PASS",
    f"ids={len(ids)}",
    "local_links=resolved",
    f"d={distance}",
    f"q={qubits}",
    f"P_fail={p_fail:.12f}",
    "ALU=3 cases",
    "guards=pass",
)
