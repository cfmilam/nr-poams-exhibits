#!/usr/bin/env python3
"""Molecular holdout for the frozen native one-eye primitive tensor.

Native one-eye action channels are translated to deterministic Gaussian
quadrature contractions solely so a mature multicentre integral engine can
evaluate shifted-centre postings.  No named basis set or molecular benchmark
enters the contractions.  Conventional names below are quarantined numerical
translation.
"""
from __future__ import annotations

import argparse
import functools
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import qr
from scipy.optimize import minimize_scalar

try:
    from pyscf import gto, scf
except ImportError as exc:
    raise SystemExit("run with the frozen PySCF-2.14 environment") from exc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import molecular_native_carbon as native
import molecular_refiling_completeness as refiling

D9 = refiling.D9
GAUSSIAN_CANDIDATES = 2.0 ** (np.arange(-32, 49, dtype=float) / 4.0)
GAUSSIAN_COUNT = 40
FIT_X, FIT_W = np.polynomial.legendre.leggauss(900)
FIT_R = (FIT_X + 1.0) * 25.0
FIT_WEIGHT = FIT_W * 25.0 * FIT_R**2


def all_roots(S, H):
    X = native.canonicalizer(S)
    values, vectors = np.linalg.eigh(X.T @ H @ X)
    return values, X @ vectors


def stationary_focks(book, result, p_count):
    Ps = result["Cs"] @ result["Cs"].T
    p = result["Cp"][:, None]
    Dp = p @ p.T
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
    pair = Jxy-Kxy if p_count == 2 else .5*Jsame+2.5*Jxy-1.5*Kxy
    Fp = book.Hp + 2*Js - Ks + pair
    return Fs, Fp


def radial_values(radials, coefficients):
    values = np.zeros_like(FIT_R)
    for coefficient, item in zip(coefficients, radials):
        polynomial = sum(v * FIT_R**k for k, v in enumerate(item.poly))
        values += coefficient * polynomial * np.exp(-item.alpha * FIT_R)
    return values


def fit_columns(radials, columns, l):
    candidates = np.stack([
        gto.gto_norm(l, exponent) * FIT_R**l
        * np.exp(-exponent * FIT_R**2)
        for exponent in GAUSSIAN_CANDIDATES
    ], axis=1)
    weight_root = np.sqrt(FIT_WEIGHT)
    weighted_candidates = candidates * weight_root[:, None]
    _, _, pivots = qr(weighted_candidates, mode="economic", pivoting=True)
    selected = np.sort(pivots[:GAUSSIAN_COUNT])
    exponents = GAUSSIAN_CANDIDATES[selected]
    primitives = candidates[:, selected]
    weighted = weighted_candidates[:, selected]
    fitted, errors = [], []
    for column in columns.T:
        target = radial_values(radials, column)
        coefficients = np.linalg.lstsq(
            weighted, target * weight_root, rcond=1e-12
        )[0]
        fit = primitives @ coefficients
        error = math.sqrt(
            np.sum(FIT_WEIGHT * (fit-target)**2)
            / np.sum(FIT_WEIGHT * target**2)
        )
        fitted.append(coefficients)
        errors.append(error)
    return exponents, np.array(fitted).T, errors


def native_channels(entry, rung):
    scalar, directed = refiling.ladders(entry)
    Cs = scalar[min(rung, len(scalar) - 1)]
    Cp = directed[min(rung, len(directed) - 1)]
    return entry.scalar_shapes, entry.directed_shapes, Cs, Cp, entry.energy


def generated_basis(entry, rung):
    s, p, Cs, Cp, native_energy = native_channels(entry, rung)
    s_exponents, gs, s_errors = fit_columns(s, Cs, 0)
    p_exponents, gp, p_errors = fit_columns(p, Cp, 1)
    shell_s = [0] + [
        [float(exponent), *map(float, row)]
        for exponent, row in zip(s_exponents, gs)
    ]
    shells = [shell_s]
    if Cp.shape[1]:
        shell_p = [1] + [
            [float(exponent), *map(float, row)]
            for exponent, row in zip(p_exponents, gp)
        ]
        shells.append(shell_p)
    return shells, native_energy, max(s_errors + p_errors)


def make_bases(strict_occupied=False):
    entries = [refiling.atom(symbol) for symbol in ("H", "C", "O")]
    rung = 0 if strict_occupied else refiling.validated_rung(entries)
    print(f"Gram-block action closure rung q={rung}")
    bases, controls = {}, {}
    for entry in entries:
        basis, native_energy, fit_error = generated_basis(entry, rung)
        bases[entry.symbol] = basis
        controls[entry.symbol] = (native_energy, fit_error)
        s_count = len(basis[0][1]) - 1
        p_count = 0 if len(basis) == 1 else len(basis[1][1]) - 1
        print(
            f"{entry.symbol} generated channels s={s_count} p={p_count} "
            f"largest radial fit error={100*fit_error:.4f}%"
        )
    return bases, controls


def total_energy(atoms, basis):
    molecule = gto.M(
        atom=atoms, basis=basis, unit="Angstrom", spin=0, charge=0, verbose=0
    )
    overlap_roots = np.linalg.eigvalsh(molecule.intor("int1e_ovlp"))
    if overlap_roots[0] < 1e-8:
        raise RuntimeError(f"molecular Gram root too small: {overlap_roots[0]:.3e}")
    mean = scf.RHF(molecule)
    mean.conv_tol = 2e-10
    mean.max_cycle = 120
    energy = mean.kernel()
    if not mean.converged:
        raise RuntimeError("stationary molecular book did not converge")
    return float(energy)


def methane_atoms(r):
    a = r / math.sqrt(3.0)
    vertices = ((1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1))
    return [("C", (0,0,0))] + [
        ("H", (a*x, a*y, a*z)) for x,y,z in vertices
    ]


def co2_atoms(r):
    return [("O", (0,0,-r)), ("C", (0,0,0)), ("O", (0,0,r))]



def ethyne_atoms(cc, ch):
    return [
        ("H", (-cc/2-ch,0,0)), ("C", (-cc/2,0,0)),
        ("C", (cc/2,0,0)), ("H", (cc/2+ch,0,0)),
    ]


def ethene_atoms(cc, ch):
    h = math.sqrt(3.0) * ch / 2.0
    return [
        ("C", (-cc/2,0,0)), ("C", (cc/2,0,0)),
        ("H", (-cc/2-ch/2,h,0)), ("H", (-cc/2-ch/2,-h,0)),
        ("H", (cc/2+ch/2,h,0)), ("H", (cc/2+ch/2,-h,0)),
    ]


def ethane_atoms(cc, ch):
    atoms = [("C", (-cc/2,0,0)), ("C", (cc/2,0,0))]
    transverse = math.sqrt(8.0/9.0)
    for side, phase in ((-1,0.0), (1,math.pi/3)):
        centre = np.array((side*cc/2,0.0,0.0))
        axial = side/3.0
        for k in range(3):
            angle = phase + 2*math.pi*k/3
            direction = np.array((
                axial, transverse*math.cos(angle), transverse*math.sin(angle)
            ))
            atoms.append(("H", tuple(centre + ch*direction)))
    return atoms


def benzene_atoms(cc, ch):
    atoms = []
    for k in range(6):
        angle = 2*math.pi*k/6
        direction = np.array((math.cos(angle),math.sin(angle),0.0))
        atoms.append(("C", tuple(cc*direction)))
        atoms.append(("H", tuple((cc+ch)*direction)))
    return atoms


def optimize_pair(name, builder, bounds, records, basis):
    cache = {}
    def ledger(cc, ch):
        key = (round(float(cc),7), round(float(ch),7))
        if key not in cache:
            cache[key] = total_energy(builder(*key), basis)
        return cache[key]
    cc = sum(bounds[0])/2
    ch = sum(bounds[1])/2
    calls = 0
    for _ in range(2):
        result = minimize_scalar(
            lambda x: ledger(x,ch), bounds=bounds[0], method="bounded",
            options={"xatol":1.5e-4,"maxiter":24},
        )
        cc = float(result.x); calls += result.nfev
        result = minimize_scalar(
            lambda x: ledger(cc,x), bounds=bounds[1], method="bounded",
            options={"xatol":1.5e-4,"maxiter":24},
        )
        ch = float(result.x); calls += result.nfev
    errors = (100*(cc/records[0]-1), 100*(ch/records[1]-1))
    print(
        f"{name} stationary lengths={cc:.8f}/{ch:.8f} A "
        f"errors={errors[0]:+.4f}%/{errors[1]:+.4f}% calls={calls}"
    )
    return (cc,ch), errors
def optimize_length(name, atom_builder, bounds, observed, basis):
    @functools.lru_cache(maxsize=None)
    def ledger(radius):
        return total_energy(atom_builder(radius), basis)
    result = minimize_scalar(
        lambda x: ledger(round(float(x), 7)),
        bounds=bounds, method="bounded",
        options={"xatol": 8e-5, "maxiter": 32},
    )
    if not result.success:
        raise RuntimeError(result.message)
    error = 100.0 * (result.x / observed - 1.0)
    print(
        f"{name} stationary length={result.x:.8f} A "
        f"record={observed:.6f} A error={error:+.4f}% calls={result.nfev}"
    )
    return float(result.x), error


def run(selected, strict_occupied=False):
    basis, _ = make_bases(strict_occupied)
    errors = []
    if selected in ("holdouts", "all"):
        methane = optimize_length(
            "methane C-H", methane_atoms, (.85, 1.35), 1.087, basis
        )
        co2 = optimize_length(
            "CO2 C-O", co2_atoms, (.95, 1.40), 1.162, basis
        )
        errors.extend((methane[1], co2[1]))
    cases = {
        "ethane": (ethane_atoms, ((1.2,1.8),(.9,1.25)), (1.535,1.094)),
        "ethyne": (ethyne_atoms, ((1.0,1.4),(.9,1.2)), (1.203,1.063)),
        "ethene": (ethene_atoms, ((1.15,1.55),(.9,1.2)), (1.339,1.086)),
        "benzene": (benzene_atoms, ((1.2,1.6),(.9,1.2)), (1.397,1.084)),
    }
    for name, (builder, bounds, records) in cases.items():
        if selected in (name, "all"):
            _, case_errors = optimize_pair(
                name, builder, bounds, records, basis
            )
            errors.extend(case_errors)
    passed = max(map(abs, errors)) < 5.0
    print("MOLECULAR BOOK", "PASS" if passed else "FAIL")
    return passed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--molecule", choices=("holdouts","ethane","ethyne","ethene","benzene","all"), default="holdouts"
    )
    parser.add_argument(
        "--strict-occupied", action="store_true",
        help="use only the frozen occupied one-eye channels for Gate 8",
    )
    args = parser.parse_args()
    raise SystemExit(
        0 if run(args.molecule, args.strict_occupied) else 2
    )


if __name__ == "__main__":
    main()
