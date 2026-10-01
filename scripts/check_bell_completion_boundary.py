#!/usr/bin/env python3
"""Finite checks for the Bell completion boundary used by the exhibit."""
from itertools import product
from math import sqrt

checks = 0
values = (-1, 1)
seen = set()
for a0, a1, b0, b1 in product(values, repeat=4):
    s = a0*b0 + a0*b1 + a1*b0 - a1*b1
    assert abs(s) <= 2
    seen.add(s)
    checks += 1
assert seen == {-2, 2}
checks += 1

strategies = list(product(values, repeat=4))
for i, left in enumerate(strategies):
    for j, right in enumerate(strategies):
        w = ((7*i + 11*j) % 101) / 100
        def score(t):
            a0, a1, b0, b1 = t
            return a0*b0 + a0*b1 + a1*b0 - a1*b1
        mixed = w*score(left) + (1-w)*score(right)
        assert abs(mixed) <= 2 + 1e-12
        checks += 1

pr = {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1): -1}
assert pr[0,0] + pr[0,1] + pr[1,0] - pr[1,1] == 4
for x, y in product((0, 1), repeat=2):
    joint = {(a, b): (0.5 if (a ^ b) == (x & y) else 0.0)
             for a, b in product((0, 1), repeat=2)}
    assert all(abs(sum(joint[a,b] for b in (0,1)) - 0.5) < 1e-12 for a in (0,1))
    assert all(abs(sum(joint[a,b] for a in (0,1)) - 0.5) < 1e-12 for b in (0,1))
    checks += 4

assert 2 < 2*sqrt(2) < 4
checks += 1
print(f"PASS: {checks} Bell-bound, convexity, PR-box, and ceiling checks")
