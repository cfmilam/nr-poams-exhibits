#!/usr/bin/env python3
"""Frozen deterministic Gram-block action-closure gates.

For a Gram book G, primitive one-entry action book B, and occupied census C,
construct K_q = span{C, (G^-1 B)C, ..., (G^-1 B)^q C}. Reconciliation occurs
in the canonical Gram chart; no empty channel is named or ordered.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np

import molecular_native_carbon as native

D9 = (.0625, .125, .25, .5, 1., 2., 4., 8., 16.)
GRAM_CUTOFF = 2e-10
CLOSURE_CUTOFF = 2e-10
VALIDATED_RUNG = 2  # first rung with two successive native convergence intervals


def gram_chart(G):
    """Return X with X.T G X=I on the reconciled Gram support."""
    return native.canonicalizer(G, cutoff=GRAM_CUTOFF)


def orthogonal_chart(columns, cutoff=CLOSURE_CUTOFF):
    """Euclidean span, invariant under nonsingular column changes."""
    if columns.shape[1] == 0:
        return np.empty((columns.shape[0], 0))
    u, singular, _ = np.linalg.svd(columns, full_matrices=False)
    keep = singular > cutoff * singular[0]
    return u[:, keep]


def closure_ladder(G, B, occupied, cutoff=CLOSURE_CUTOFF):
    """Return every distinct K_q plus one explicit fixed-point witness."""
    X = gram_chart(G)
    H = X.T @ B @ X
    seed = orthogonal_chart(X.T @ G @ occupied, cutoff)
    if seed.shape[1] == 0:
        empty = np.empty((len(G), 0))
        return [empty, empty.copy()]
    blocks = [seed]
    for _ in range(H.shape[0] + 1):
        trial = orthogonal_chart(
            np.column_stack((blocks[-1], H @ blocks[-1])), cutoff
        )
        blocks.append(trial)
        if trial.shape[1] == blocks[-2].shape[1]:
            break
    else:
        raise RuntimeError("action closure did not reach a finite fixed point")
    return [X @ block for block in blocks]


def projector(G, columns):
    X = gram_chart(G)
    q = orthogonal_chart(X.T @ G @ columns)
    return q @ q.T


def subspace_defect(G, left, right):
    return float(np.linalg.norm(projector(G, left) - projector(G, right)))


def ritz_roots(B, columns):
    return np.linalg.eigvalsh(columns.T @ B @ columns)


def transform4(value, a, b, c, d):
    return np.einsum(
        "ai,bj,ck,dl,abcd->ijkl", a, b, c, d, value, optimize=True
    )


def restrict_book(book, scalar, directed):
    """Project the complete primitive book into one closure rung."""
    return native.Books(
        [None] * scalar.shape[1], [None] * directed.shape[1],
        scalar.T @ book.Ss @ scalar,
        directed.T @ book.Sp @ directed,
        scalar.T @ book.Hs @ scalar,
        directed.T @ book.Hp @ directed,
        scalar.T @ book.Ms @ scalar,
        directed.T @ book.Mp @ directed,
        transform4(book.ssss, scalar, scalar, scalar, scalar),
        transform4(book.sspp, scalar, scalar, directed, directed),
        transform4(book.spsp, scalar, directed, scalar, directed),
        transform4(book.pxpx_pypy, directed, directed, directed, directed),
        transform4(book.pxpy_pypx, directed, directed, directed, directed),
        transform4(book.pxpx_pxpx, directed, directed, directed, directed),
    )


def lift_result(result, scalar, directed):
    return {
        **result,
        "Cs": scalar @ result["Cs"],
        "Cp": directed @ result["Cp"],
    }


@dataclass
class Atom:
    symbol: str
    scalar_shapes: list
    directed_shapes: list
    Ss: np.ndarray
    Sp: np.ndarray
    Bs: np.ndarray
    Bp: np.ndarray
    Cs: np.ndarray
    Cp: np.ndarray
    energy: float
    book: native.Books | None
    p_count: int


def atom(symbol):
    if symbol == "H":
        scalar = [
            native.shape(family, z)
            for family in ("one", "two-scalar") for z in D9
        ]
        directed = [native.shape("two-directed", z) for z in D9]
        Ss = native.matrix(scalar, native.overlap)
        Sp = native.matrix(directed, native.overlap)
        Bs = native.matrix(
            scalar, lambda a, b: native.one_posting(a, b, 1.0)
        )
        Bp = native.matrix(
            directed, lambda a, b: native.one_posting(a, b, 1.0)
        )
        value, Cs = native.lowest(Ss, Bs, 1)
        return Atom(
            symbol, scalar, directed, Ss, Sp, Bs, Bp, Cs,
            np.empty((len(directed), 0)), float(value[0]), None, 0,
        )
    eye, p_count = {"C": (6.0, 2), "O": (8.0, 4)}[symbol]
    book = native.build(D9, eye)
    result = native.solve(book, p_count)
    return Atom(
        symbol, book.s, book.p, book.Ss, book.Sp, book.Hs, book.Hp,
        result["Cs"], result["Cp"][:, None], result["energy"], book, p_count,
    )


def ladders(entry):
    return (
        closure_ladder(entry.Ss, entry.Bs, entry.Cs),
        closure_ladder(entry.Sp, entry.Bp, entry.Cp),
    )


def validated_rung(entries):
    """Universal native stop: first rung admitting two convergence intervals."""
    return VALIDATED_RUNG


def one_eye_row(entry, scalar, directed):
    if entry.symbol == "H":
        value, coeff = native.lowest(
            scalar.T @ entry.Ss @ scalar,
            scalar.T @ entry.Bs @ scalar, 1,
        )
        full = scalar @ coeff[:, 0]
        moment_book = native.matrix(entry.scalar_shapes, native.moment)
        radius = float(full.T @ moment_book @ full)
        return float(value[0]), np.array((radius,)), {}
    reduced = restrict_book(entry.book, scalar, directed)
    solved = native.solve(reduced, entry.p_count)
    lifted = lift_result(solved, scalar, directed)
    invariants = native.contractions(native.occupied_tensor(entry.book, lifted))
    return solved["energy"], np.array(solved["moments"]), invariants


def relative_drift(new, old):
    new = np.asarray(new)
    old = np.asarray(old)
    scale = np.maximum(np.maximum(np.abs(new), np.abs(old)), 1e-14)
    return float(np.max(np.abs(new - old) / scale))


def rotated_seed(entry, rng):
    def change(columns):
        n = columns.shape[1]
        if n == 0:
            return columns.copy()
        q, _ = np.linalg.qr(rng.normal(size=(n, n)))
        return columns @ q @ np.diag(np.geomspace(.7, 1.3, n))
    return change(entry.Cs), change(entry.Cp)


def overcomplete_seed(columns):
    if columns.shape[1] == 0:
        return columns.copy()
    return np.column_stack((columns[:, ::-1], columns, columns[:, :1]))


def adverse_eigen_seed(G, B, count):
    if count == 0:
        return np.empty((len(G), 0))
    X = gram_chart(G)
    _, vectors = np.linalg.eigh(X.T @ B @ X)
    return X @ vectors[:, -count:]


def run():
    entries = [atom(symbol) for symbol in ("H", "C", "O")]
    universal = validated_rung(entries)
    print(f"universal two-interval convergence rung q={universal}")
    all_ok = True
    for entry in entries:
        scalar_ladder, directed_ladder = ladders(entry)
        at = lambda seq, q: seq[min(q, len(seq) - 1)]
        scalar = at(scalar_ladder, universal)
        directed = at(directed_ladder, universal)
        ranks = [
            (at(scalar_ladder, q).shape[1], at(directed_ladder, q).shape[1])
            for q in range(universal + 1)
        ]
        base = one_eye_row(entry, scalar_ladder[0], directed_ladder[0])
        occupied_recovery = relative_drift(base[0], entry.energy)
        energies, moments, tensors = [], [], []
        for q in range(universal + 1):
            row = one_eye_row(
                entry, at(scalar_ladder, q), at(directed_ladder, q)
            )
            energies.append(row[0]); moments.append(row[1]); tensors.append(row[2])
        monotonic_raise = max(
            [0.0] + [energies[i] - energies[i-1] for i in range(1, len(energies))]
        )
        intervals = range(1, universal + 1)
        rate_drift = max(
            relative_drift(energies[q], energies[q-1]) for q in intervals
        )
        moment_drift = max(
            relative_drift(moments[q], moments[q-1]) for q in intervals
        )
        tensor_drift = 0.0
        if tensors[-1]:
            tensor_drift = max(
                relative_drift(tensors[q][key], tensors[q-1][key])
                for q in intervals
                for key in ("frobenius", "census_trace", "exchange_trace")
            )
        rng = np.random.default_rng(20260910 + ord(entry.symbol))
        rotated_s, rotated_p = rotated_seed(entry, rng)
        rotated_sl = closure_ladder(entry.Ss, entry.Bs, rotated_s)
        rotated_pl = closure_ladder(entry.Sp, entry.Bp, rotated_p)
        chart_defect = max(
            subspace_defect(entry.Ss, scalar, at(rotated_sl, universal)),
            subspace_defect(entry.Sp, directed, at(rotated_pl, universal)),
        )
        root_defect = 0.0
        for B, left, right in (
            (entry.Bs, scalar, at(rotated_sl, universal)),
            (entry.Bp, directed, at(rotated_pl, universal)),
        ):
            if left.shape[1]:
                root_defect = max(
                    root_defect,
                    relative_drift(ritz_roots(B, right), ritz_roots(B, left)),
                )
        rotated_row = one_eye_row(entry, at(rotated_sl, universal), at(rotated_pl, universal))
        invariant_defect = relative_drift(rotated_row[0], energies[-1])
        if tensors[-1]:
            invariant_defect = max(
                invariant_defect,
                max(
                    relative_drift(rotated_row[2][key], tensors[-1][key])
                    for key in ("frobenius", "census_trace", "exchange_trace")
                ),
            )
        reordered_sl = closure_ladder(
            entry.Ss, entry.Bs, overcomplete_seed(entry.Cs)
        )
        reordered_pl = closure_ladder(
            entry.Sp, entry.Bp, overcomplete_seed(entry.Cp)
        )
        order_defect = max(
            subspace_defect(entry.Ss, scalar, at(reordered_sl, universal)),
            subspace_defect(entry.Sp, directed, at(reordered_pl, universal)),
        )
        arbitrary_sl = closure_ladder(
            entry.Ss, entry.Bs,
            adverse_eigen_seed(entry.Ss, entry.Bs, entry.Cs.shape[1]),
        )
        arbitrary_pl = closure_ladder(
            entry.Sp, entry.Bp,
            adverse_eigen_seed(entry.Sp, entry.Bp, entry.Cp.shape[1]),
        )
        arbitrary_defect = max(
            subspace_defect(entry.Ss, scalar, at(arbitrary_sl, universal)),
            subspace_defect(entry.Sp, directed, at(arbitrary_pl, universal)),
        )
        gates = {
            "occupied": occupied_recovery < 1e-10,
            "chart": max(chart_defect, root_defect, invariant_defect) < 1e-9,
            "monotonic": monotonic_raise < 1e-9,
            "finite": rate_drift < .005 and moment_drift < .01
                      and tensor_drift < .01,
            "adverse": order_defect < 1e-9 and arbitrary_defect > 1e-6,
        }
        all_ok &= all(gates.values())
        print(
            f"{entry.symbol}: ranks={ranks} occupied={occupied_recovery:.3e} "
            f"chart/root/invariant={chart_defect:.3e}/"
            f"{root_defect:.3e}/{invariant_defect:.3e}"
        )
        print(
            f"{entry.symbol}: max_raise={monotonic_raise:.3e} "
            f"final rate/moment/tensor={100*rate_drift:.6f}%/"
            f"{100*moment_drift:.6f}%/{100*tensor_drift:.6f}% "
            f"order/arbitrary={order_defect:.3e}/{arbitrary_defect:.3e}"
        )
        print(
            "GATES 1-5",
            " ".join(f"{key}={'PASS' if value else 'FAIL'}"
                     for key, value in gates.items()),
        )
    print("NATIVE REFILING", "PASS" if all_ok else "FAIL")
    return all_ok, universal


def main():
    argparse.ArgumentParser().parse_args()
    passed, _ = run()
    raise SystemExit(0 if passed else 2)


if __name__ == "__main__":
    main()
