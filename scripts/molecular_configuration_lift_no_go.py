#!/usr/bin/env python3
"""Test whether the presently stated one-entry axioms uniquely force a many-entry book.

This is a structural counterexample, not a molecular model.  It uses three admissible
orientation-and-sense slots and the three two-entry configurations 01, 02, 12.
Duplicate occupancy is absent by construction (the exterior census).
"""
from __future__ import annotations

import itertools
import numpy as np


def exterior_two_book(g: np.ndarray, b: np.ndarray):
    """Induced no-extra-posting two-entry books in the ordered wedge chart."""
    pairs = list(itertools.combinations(range(len(g)), 2))
    n = len(pairs)
    g2 = np.empty((n, n))
    b2 = np.empty((n, n))
    for p, (i, j) in enumerate(pairs):
        for q, (k, l) in enumerate(pairs):
            g2[p, q] = g[i, k] * g[j, l] - g[i, l] * g[j, k]
            b2[p, q] = (
                b[i, k] * g[j, l] + g[i, k] * b[j, l]
                - b[i, l] * g[j, k] - g[i, l] * b[j, k]
            )
    return pairs, g2, b2


def roots(b: np.ndarray, g: np.ndarray):
    l = np.linalg.cholesky(g)
    a = np.linalg.solve(l, b)
    a = np.linalg.solve(l, a.T).T
    return np.linalg.eigvalsh((a + a.T) / 2)


if __name__ == "__main__":
    # One-entry book: fixed input shared by both lifts.
    g = np.eye(3)
    b = np.diag((-1.0, -0.4, 0.3))
    pairs, g2, base = exterior_two_book(g, b)

    # A distinct, finite pair-posting reconciliation on the same admissible
    # configurations. It introduces no duplicate slot and remains a valid symmetric
    # ledger. Existing axioms do not specify or forbid its coefficient.
    pair_posting = np.diag((0.00, 0.25, -0.15))
    alt = base + pair_posting

    r0 = roots(base, g2)
    r1 = roots(alt, g2)
    print("configurations:", pairs)
    print("induced roots:", np.array2string(r0, precision=9))
    print("alternate roots:", np.array2string(r1, precision=9))
    print("root difference:", np.array2string(r1-r0, precision=9))
    assert np.linalg.det(g2) > 0
    assert np.max(np.abs(r1-r0)) > 0.1

    # Both books retain generalized-root invariance under an arbitrary nonsingular
    # re-description of the configuration chart.
    q = np.array(((1.0, .2, -.1), (.1, 1.1, .3), (-.2, .1, .9)))
    for book in (base, alt):
        rq = roots(q.T @ book @ q, q.T @ g2 @ q)
        assert np.max(np.abs(rq-roots(book, g2))) < 1e-12

    print("VERDICT: EXISTING MANY-ENTRY PREMISES DO NOT UNIQUELY FIX B_config")
