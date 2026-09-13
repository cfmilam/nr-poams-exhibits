#!/usr/bin/env python3
"""Finite pair-only lift from slot postings to whole occupied configurations.

The occupied configurations are ordered exterior census entries.  `h[p,q]` is a
one-entry posting and `v[p,q,r,s]` is the primitive two-entry posting.  No
many-entry coefficient appears: applying those primitive postings to every
occupied slot or occupied pair constructs the configuration book.
"""
from __future__ import annotations

import itertools
import numpy as np


def apply_remove(state: int, q: int):
    if not state >> q & 1:
        return None
    sign = -1 if (state & ((1 << q) - 1)).bit_count() & 1 else 1
    return sign, state ^ (1 << q)


def apply_add(state: int, p: int):
    if state >> p & 1:
        return None
    sign = -1 if (state & ((1 << p) - 1)).bit_count() & 1 else 1
    return sign, state | (1 << p)


def configurations(nslots: int, nentries: int):
    return [sum(1 << i for i in occ) for occ in itertools.combinations(range(nslots), nentries)]


def lift(h: np.ndarray, v: np.ndarray, nentries: int):
    """Construct B_config from primitive one- and two-entry postings."""
    n = len(h); states = configurations(n, nentries); index = {s:i for i,s in enumerate(states)}
    out = np.zeros((len(states), len(states)))
    for col, state in enumerate(states):
        for p in range(n):
            for q in range(n):
                aq = apply_remove(state, q)
                if aq is None: continue
                cp = apply_add(aq[1], p)
                if cp is not None:
                    out[index[cp[1]], col] += aq[0]*cp[0]*h[p,q]
        for p,q,r,s in itertools.product(range(n), repeat=4):
            a_s = apply_remove(state, s)
            if a_s is None: continue
            a_r = apply_remove(a_s[1], r)
            if a_r is None: continue
            c_q = apply_add(a_r[1], q)
            if c_q is None: continue
            c_p = apply_add(c_q[1], p)
            if c_p is not None:
                out[index[c_p[1]], col] += .5*a_s[0]*a_r[0]*c_q[0]*c_p[0]*v[p,q,r,s]
    return states, (out+out.T)/2


if __name__ == "__main__":
    rng=np.random.default_rng(20260910); n=5
    h=rng.normal(size=(n,n));h=(h+h.T)/2
    raw=rng.normal(size=(n,n,n,n))
    # A real two-ended pair book: interchange of the two entries and reversal
    # of bra/ket descriptions leave the posting unchanged.
    v=(raw+raw.transpose(1,0,3,2)+raw.transpose(2,3,0,1)+raw.transpose(3,2,1,0))/4
    states,b=lift(h,v,2)
    print(f"slots={n} admissible_two_entry_configurations={len(states)}")
    print(f"configuration_book_symmetry={np.max(np.abs(b-b.T)):.3e}")
    assert len(states)==10 and np.max(np.abs(b-b.T))<1e-12

    # Duplicate occupation is not an admissible state: C(5,2), never 5^2.
    assert all(s.bit_count()==2 for s in states)

    # With no pair posting, diagonal roots are exactly sums of occupied
    # one-entry postings. This is the required separated-entry limit.
    diag=np.diag(np.arange(n,dtype=float));_,sep=lift(diag,np.zeros_like(v),2)
    expected=sorted(i+j for i,j in itertools.combinations(range(n),2))
    measured=sorted(np.linalg.eigvalsh(sep))
    print("separated-limit roots:", measured)
    assert np.max(np.abs(np.array(expected)-measured))<1e-12
    print("ALL PRIMITIVE-PAIR LIFT GATES PASS")
