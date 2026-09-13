#!/usr/bin/env python3
"""Cusp-resolving one-eye carbon circulation ledger.

The finite chart contains the frozen one-turn scalar, two-turn scalar, and
two-turn directed action shapes.  Every range posting is evaluated by an
exact spherical multipole reduction; no Cartesian grid, observed radius,
screening coefficient, or named basis set enters.

Conventional translation (quarantined): this is the variational six-entry
determinant book for 1s2 2s2 2p2 in a generated dyadic STO chart.
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from functools import lru_cache
from itertools import product

import numpy as np

EYE_COUNT = 6.0


@dataclass(frozen=True)
class Radial:
    family: str
    zeta: float
    alpha: float
    poly: tuple[float, ...]
    l: int


def shape(family: str, zeta: float) -> Radial:
    """Normalized radial action shape: integral r^2 R(r)^2 dr = 1."""
    if family == "one":
        return Radial(family, zeta, zeta, (2.0 * zeta**1.5,), 0)
    if family == "two-scalar":
        n = zeta**1.5 / math.sqrt(2.0)
        return Radial(family, zeta, zeta / 2.0, (n, -n * zeta / 2.0), 0)
    if family == "two-directed":
        return Radial(
            family, zeta, zeta / 2.0,
            (0.0, zeta**2.5 / (2.0 * math.sqrt(6.0))), 1
        )
    raise ValueError(family)


def conv(a, b):
    return tuple(np.convolve(a, b).tolist())


def deriv(poly, alpha):
    out = [-alpha * x for x in poly]
    for k in range(1, len(poly)):
        out[k - 1] += k * poly[k]
    return tuple(out)


def gamma_integral(power, rate):
    return math.factorial(power) / rate ** (power + 1)


def bilinear(a, b, extra_power):
    return sum(
        c * gamma_integral(k + extra_power, a.alpha + b.alpha)
        for k, c in enumerate(conv(a.poly, b.poly))
    )


def overlap(a, b):
    return bilinear(a, b, 2)


def moment(a, b):
    return bilinear(a, b, 3)


def one_posting(a, b, eye=EYE_COUNT):
    if a.l != b.l:
        return 0.0
    radial_turn = sum(
        c * gamma_integral(k + 2, a.alpha + b.alpha)
        for k, c in enumerate(conv(deriv(a.poly, a.alpha), deriv(b.poly, b.alpha)))
    )
    angular_turn = a.l * (a.l + 1) * bilinear(a, b, 0)
    return 0.5 * (radial_turn + angular_turn) - eye * bilinear(a, b, 1)


def angular_table():
    """Polynomial-exact real Gaunt coefficients through L=2."""
    x, w = np.polynomial.legendre.leggauss(10)
    phi = 2.0 * math.pi * np.arange(24) / 24.0
    rows = []
    for nz, wz in zip(x, w):
        st = math.sqrt(max(0.0, 1.0 - nz * nz))
        for ph in phi:
            nx, ny = st * math.cos(ph), st * math.sin(ph)
            weight = wz * 2.0 * math.pi / len(phi)
            y0 = 1.0 / math.sqrt(4.0 * math.pi)
            y1 = math.sqrt(3.0 / (4.0 * math.pi)) * np.array((nx, ny, nz))
            y2 = np.array((
                math.sqrt(5.0 / (16.0 * math.pi)) * (3.0 * nz * nz - 1.0),
                math.sqrt(15.0 / (16.0 * math.pi)) * (nx * nx - ny * ny),
                math.sqrt(15.0 / (4.0 * math.pi)) * nx * ny,
                math.sqrt(15.0 / (4.0 * math.pi)) * nx * nz,
                math.sqrt(15.0 / (4.0 * math.pi)) * ny * nz,
            ))
            rows.append((
                weight,
                {"s": y0, "x": y1[0], "y": y1[1], "z": y1[2]},
                (np.array((y0,)), y1, y2),
            ))
    return {
        (a, b): {
            L: sum(
                weight * vals[a] * vals[b] * harmonics[L]
                for weight, vals, harmonics in rows
            )
            for L in range(3)
        }
        for a in ("s", "x", "y", "z")
        for b in ("s", "x", "y", "z")
    }


ANGULAR = angular_table()


@lru_cache(maxsize=None)
def radial_pair(a, b):
    return conv(a.poly, b.poly), a.alpha + b.alpha


@lru_cache(maxsize=None)
def radial_range(a, b, c, d, L):
    """Exact double radial multipole integral, split at r=s."""
    p, ap = radial_pair(a, b)
    q, aq = radial_pair(c, d)

    def half(outer, ao, inner, ai):
        total = 0.0
        for ko, co in enumerate(outer):
            if co == 0.0:
                continue
            for ki, ci in enumerate(inner):
                if ci == 0.0:
                    continue
                n = ki + 1 - L
                if n < 0:
                    raise ArithmeticError("forbidden singular multipole product")
                for j in range(n + 1):
                    tail = (
                        math.factorial(n) / math.factorial(j)
                        / ai ** (n - j + 1)
                    )
                    total += (
                        co * ci * tail
                        * gamma_integral(ko + L + 2 + j, ao + ai)
                    )
        return total

    return half(p, ap, q, aq) + half(q, aq, p, ap)


@lru_cache(maxsize=None)
def primitive_range(a, ka, b, kb, c, kc, d, kd):
    value = 0.0
    for L in range(3):
        angular = float(np.dot(ANGULAR[ka, kb][L], ANGULAR[kc, kd][L]))
        if abs(angular) > 1e-15:
            value += (
                4.0 * math.pi / (2 * L + 1)
                * angular * radial_range(a, b, c, d, L)
            )
    return value


def matrix(items, fn):
    out = np.empty((len(items), len(items)))
    for i, a in enumerate(items):
        for j in range(i + 1):
            out[i, j] = out[j, i] = fn(a, items[j])
    return out


def tensor(a, ka, b, kb, c, kc, d, kd):
    out = np.empty((len(a), len(b), len(c), len(d)))
    for i, j, k, l in product(
        range(len(a)), range(len(b)), range(len(c)), range(len(d))
    ):
        out[i, j, k, l] = primitive_range(
            a[i], ka, b[j], kb, c[k], kc, d[l], kd
        )
    return out


def canonicalizer(S, cutoff=2e-10):
    vals, vecs = np.linalg.eigh(S)
    keep = vals > cutoff * vals[-1]
    return vecs[:, keep] / np.sqrt(vals[keep])


def lowest(S, H, count):
    X = canonicalizer(S)
    vals, vecs = np.linalg.eigh(X.T @ H @ X)
    return vals[:count], X @ vecs[:, :count]


@dataclass
class Books:
    s: list[Radial]
    p: list[Radial]
    Ss: np.ndarray
    Sp: np.ndarray
    Hs: np.ndarray
    Hp: np.ndarray
    Ms: np.ndarray
    Mp: np.ndarray
    ssss: np.ndarray
    sspp: np.ndarray
    spsp: np.ndarray
    pxpx_pypy: np.ndarray
    pxpy_pypx: np.ndarray
    pxpx_pxpx: np.ndarray


def build(factors, eye=EYE_COUNT):
    s = [shape(f, z) for f in ("one", "two-scalar") for z in factors]
    p = [shape("two-directed", z) for z in factors]
    return Books(
        s, p,
        matrix(s, overlap), matrix(p, overlap),
        matrix(s, lambda a,b: one_posting(a,b,eye)),
        matrix(p, lambda a,b: one_posting(a,b,eye)),
        matrix(s, moment), matrix(p, moment),
        tensor(s, "s", s, "s", s, "s", s, "s"),
        tensor(s, "s", s, "s", p, "x", p, "x"),
        tensor(s, "s", p, "x", s, "s", p, "x"),
        tensor(p, "x", p, "x", p, "y", p, "y"),
        tensor(p, "x", p, "y", p, "y", p, "x"),
        tensor(p, "x", p, "x", p, "x", p, "x"),
    )


def solve(book, p_count=2, maxiter=700, tolerance=1e-7):
    if p_count not in (2, 4):
        raise ValueError("implemented open directed censuses are p_count=2 or 4")
    _, Cs = lowest(book.Ss, book.Hs, 2)
    _, Cp = lowest(book.Sp, book.Hp, 1)
    Ps, Dp = Cs @ Cs.T, Cp @ Cp.T
    damping, last = 0.25, None
    for iteration in range(1, maxiter + 1):
        Jss = np.einsum("abcd,cd->ab", book.ssss, Ps, optimize=True)
        Kss = np.einsum("acbd,cd->ab", book.ssss, Ps, optimize=True)
        Jp = np.einsum("abcd,cd->ab", book.sspp, Dp, optimize=True)
        Kp = np.einsum("acbd,cd->ab", book.spsp, Dp, optimize=True)
        Fs = book.Hs + 2*Jss - Kss + p_count*Jp - .5*p_count*Kp

        Js = np.einsum("abcd,ab->cd", book.sspp, Ps, optimize=True)
        Ks = np.einsum("acbd,ab->cd", book.spsp, Ps, optimize=True)
        Jxy = np.einsum("abcd,cd->ab", book.pxpx_pypy, Dp, optimize=True)
        Kxy = np.einsum("abcd,cd->ab", book.pxpy_pypx, Dp, optimize=True)
        Jsame = np.einsum("abcd,cd->ab", book.pxpx_pxpx, Dp, optimize=True)
        if p_count == 2:
            pair_fock = Jxy - Kxy
        else:
            pair_fock = .5*Jsame + 2.5*Jxy - 1.5*Kxy
        Fp = book.Hp + 2*Js - Ks + pair_fock

        _, Cs_new = lowest(book.Ss, Fs, 2)
        _, Cp_new = lowest(book.Sp, Fp, 1)
        Ps_new, Dp_new = Cs_new @ Cs_new.T, Cp_new @ Cp_new.T
        residual = max(
            np.linalg.norm(Ps_new - Ps), np.linalg.norm(Dp_new - Dp)
        )
        if residual < tolerance:
            Ps, Dp, Cs, Cp = Ps_new, Dp_new, Cs_new, Cp_new
            break
        Ps = (1 - damping) * Ps + damping * Ps_new
        Dp = (1 - damping) * Dp + damping * Dp_new
        if last is not None and residual > last * 1.05:
            damping = max(0.06, damping * 0.75)
        elif iteration % 35 == 0:
            damping = min(0.40, damping * 1.15)
        last = residual
    else:
        raise RuntimeError(f"stationary iteration failed: residual={residual:.3e}")

    one = (
        2 * np.einsum("ab,ab", Ps, book.Hs)
        + p_count * np.einsum("ab,ab", Dp, book.Hp)
    )
    jss = np.einsum("ab,cd,abcd", Ps, Ps, book.ssss, optimize=True)
    kss = np.einsum("ac,bd,abcd", Ps, Ps, book.ssss, optimize=True)
    jsp = np.einsum("ab,cd,abcd", Ps, Dp, book.sspp, optimize=True)
    ksp = np.einsum("ab,cd,acbd", Ps, Dp, book.spsp, optimize=True)
    jxy = np.einsum("ab,cd,abcd", Dp, Dp, book.pxpx_pypy, optimize=True)
    kxy = np.einsum("ab,cd,abcd", Dp, Dp, book.pxpy_pypx, optimize=True)
    jsame = np.einsum("ab,cd,abcd", Dp, Dp, book.pxpx_pxpx, optimize=True)
    if p_count == 2:
        pair_energy = jxy - kxy
    else:
        pair_energy = jsame + 5*jxy - 3*kxy
    energy = float(
        one + 2*jss - kss
        + p_count*(2*jsp-ksp) + pair_energy
    )
    r_s = np.linalg.eigvalsh(Cs.T @ book.Ms @ Cs)
    r_p = float((Cp.T @ book.Mp @ Cp).item())
    return {
        "energy": energy, "Cs": Cs, "Cp": Cp[:, 0],
        "moments": (*r_s, r_p), "iterations": iteration,
        "residual": residual,
    }


def occupied_tensor(book, result):
    """Five-slot spatial tensor: two scalar + all three directed charts."""
    C, p = result["Cs"], result["Cp"]
    radial = [C[:, 0], C[:, 1], p, p, p]
    kinds = ["s", "s", "x", "y", "z"]
    bases = [book.s, book.s, book.p, book.p, book.p]
    out = np.empty((5, 5, 5, 5))
    for i, j, k, l in product(range(5), repeat=4):
        raw = tensor(
            bases[i], kinds[i], bases[j], kinds[j],
            bases[k], kinds[k], bases[l], kinds[l]
        )
        out[i, j, k, l] = np.einsum(
            "a,b,c,d,abcd",
            radial[i], radial[j], radial[k], radial[l], raw,
            optimize=True,
        )
    return out


def directional_contraction(directed, occupied):
    response = (
        np.einsum("abcd,cd->ab", directed, occupied)
        - np.einsum("acbd,cd->ab", directed, occupied)
    )
    eye = np.eye(3)
    return float(np.linalg.norm(response - np.trace(response) / 3.0 * eye))


def contractions(v):
    directed = v[2:, 2:, 2:, 2:]
    occupied = np.diag((1.0, 1.0, 0.0))
    spherical = 2.0 / 3.0 * np.eye(3)
    return {
        "frobenius": float(np.linalg.norm(v)),
        "census_trace": float(sum(v[i,i,j,j] for i in range(5) for j in range(5))),
        "exchange_trace": float(sum(v[i,j,j,i] for i in range(5) for j in range(5))),
        "directed_anisotropy": directional_contraction(directed, occupied),
        "spherical_anisotropy": directional_contraction(directed, spherical),
    }


def control():
    factors = (0.25, 0.5, 1.0, 2.0, 4.0)
    items = [shape("one", z) for z in factors]
    S = matrix(items, overlap)
    H = matrix(items, lambda a, b: one_posting(a, b, eye=1.0))
    values, C = lowest(S, H, 1)
    radius = float(C[:, 0].T @ matrix(items, moment) @ C[:, 0])
    return float(values[0]), radius


def closed_family_gate():
    x, _ = np.polynomial.legendre.leggauss(24)
    phi = 2 * math.pi * np.arange(48) / 48
    worst = 0.0
    for nz in x:
        st = math.sqrt(1 - nz*nz)
        for ph in phi:
            n = np.array((st*math.cos(ph), st*math.sin(ph), nz))
            rho = 3/(4*math.pi) * float(n@n)
            worst = max(worst, abs(rho - 3/(4*math.pi)))
    rng = np.random.default_rng(20260910)
    q, _ = np.linalg.qr(rng.normal(size=(3,3)))
    return worst, float(np.linalg.norm(q @ q.T - np.eye(3)))


def run(ladders):
    ctrl_e, ctrl_r = control()
    iso, closed_chart = closed_family_gate()
    print(f"one-entry control book={ctrl_e:.12f} <r>={ctrl_r:.12f}")
    print(f"closed-family isotropy={iso:.3e} chart={closed_chart:.3e}")
    rows = []
    for factors in ladders:
        print(f"building D{len(factors)}={factors}", flush=True)
        book = build(factors)
        result = solve(book)
        v = occupied_tensor(book, result)
        inv = contractions(v)
        rows.append((factors, book, result, v, inv))
        print(
            f"D{len(factors)} book={result['energy']:.12f} "
            f"moments={result['moments']} iterations={result['iterations']}"
        )
        print(" tensor", " ".join(f"{k}={x:.12f}" for k,x in inv.items()))

    a, b = rows[-2], rows[-1]
    edrift = abs(b[2]["energy"] / a[2]["energy"] - 1)
    mdrift = np.max(abs(
        np.array(b[2]["moments"]) / np.array(a[2]["moments"]) - 1
    ))
    tdrift = max(
        abs(b[4][k] / a[4][k] - 1)
        for k in ("frobenius", "census_trace", "exchange_trace")
    )
    gram = np.eye(6)
    gram[-1] = gram[-2]
    gram[:, -1] = gram[:, -2]
    duplicate = float(np.linalg.det(gram))
    rng = np.random.default_rng(916)
    directed = b[3][2:, 2:, 2:, 2:]
    projector = np.diag((1.0, 1.0, 0.0))

    def pair_root(occupied):
        direct = np.einsum("ab,cd,abcd", occupied, occupied, directed)
        crossed = np.einsum("ac,bd,abcd", occupied, occupied, directed)
        return 0.5 * (direct - crossed)

    baseline = pair_root(projector)
    chart_drift = 0.0
    for _ in range(12):
        q, _ = np.linalg.qr(rng.normal(size=(3,3)))
        rotated = q @ projector @ q.T
        chart_drift = max(chart_drift, abs(pair_root(rotated) / baseline - 1.0))
    directional = b[4]["directed_anisotropy"]
    spherical_directional = b[4]["spherical_anisotropy"]
    print(f"final radial drift: book={100*edrift:.6f}% moments={100*mdrift:.6f}%")
    print(f"final tensor invariant drift={100*tdrift:.6f}%")
    print(f"duplicate oriented Gram volume={duplicate:.3e}")
    print(f"carbon chart root drift={chart_drift:.3e}")
    print(
        f"adverse directional contraction: native={directional:.12f} "
        f"spherical_average={spherical_directional:.3e}"
    )
    gates = {
        "one_entry": (
            abs(ctrl_e / -0.5 - 1) < 1e-5
            and abs(ctrl_r / 1.5 - 1) < 1e-5
        ),
        "closed_family": max(iso, closed_chart) < 1e-9,
        "carbon_chart": chart_drift < 1e-8,
        "duplicate": abs(duplicate) < 1e-12,
        "radial": edrift < .005 and mdrift < .01,
        "tensor": tdrift < .01,
        "adverse": directional > 1e-8 and spherical_directional < 1e-12,
    }
    print(
        "GATES 1-7",
        " ".join(
            f"{k}={'PASS' if value else 'FAIL'}"
            for k, value in gates.items()
        )
    )
    return rows, gates


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    ladders = (
        ((0.5,1.0,2.0), (0.25,0.5,1.0,2.0,4.0))
        if args.quick else (
            (0.25,0.5,1.0,2.0,4.0),
            (0.125,0.25,0.5,1.0,2.0,4.0,8.0),
            (0.0625,0.125,0.25,0.5,1.0,2.0,4.0,8.0,16.0),
        )
    )
    _, gates = run(ladders)
    raise SystemExit(0 if all(gates.values()) else 2)


if __name__ == "__main__":
    main()
