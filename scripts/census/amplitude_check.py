#!/usr/bin/env python3
"""
Matched-bias quantification (2026-08-28) — no new census runs needed.
The 08-27 verdict: class D2 rank tracks clustering amplitude (b·D(z)) EXCEPT QSO,
which reads ~+0.2 above its station. Make that quantitative INTERNALLY (no literature b):

  rho_X(r) = Gamma*_X(r) / nbar_X        (conditional overdensity around galaxies)
  rho - 1  ~ (b·D)^2 · xibar_m(r)        (extant reading of the amplitude)

nbar_X = N_used / V(window); V from the class's own randoms: footprint Omega via the
pipeline's 0.5-deg cell method + lcdm comoving shell volume. Same map (lcdm) as the D2s.
Output: per class, rho(r)-1 at pivots, implied (bD)/(bD)_LRG, measured D2 -> the
amplitude station vs dimension table, and QSO's residual off the relation.
"""
import json, numpy as np
from astropy.io import fits
from census_run import d_lcdm, load_cat

C = 0.5  # deg cell, matches estimator

def omega_sr(ra, dec):
    cell = np.radians(C)
    dec_r, ra_r = np.radians(dec), np.radians(ra)
    iy = np.floor(dec_r/cell).astype(np.int64)
    ix = np.floor(ra_r*np.cos(dec_r)/cell).astype(np.int64)
    keys = ix*100000 + iy
    uniq, ucnt = np.unique(keys, return_counts=True)
    return (ucnt >= 5).sum()*cell*cell

CLASSES = [
    # (label, json, dat fits, ran fits, zmin, zmax)
    ("BGS",  "results/BGS_BRIGHT-21.5_NGC_census.json", "data/BGS_BRIGHT-21.5_NGC_clustering.dat.fits", "data/BGS_BRIGHT-21.5_NGC_0_clustering.ran.fits", 0.05, 0.40),
    ("LRG",  "results/LRG_NGC_z08-11_census.json",      "data/LRG_NGC_clustering.dat.fits",             "data/LRG_NGC_0_clustering.ran.fits",             0.8,  1.1),
    ("QSO",  "results/QSO_NGC_z08-11_census.json",      "data/QSO_NGC_clustering.dat.fits",             "data/QSO_NGC_0_clustering.ran.fits",             0.8,  1.1),
    ("ELG",  "results/ELG_LOPnotqso_NGC_z08-11_census.json", "data/ELG_LOPnotqso_NGC_clustering.dat.fits", "data/ELG_LOPnotqso_NGC_0_clustering.ran.fits", 0.8, 1.1),
]
PIVOTS = [20.0, 30.0, 45.0, 70.0]

rows = []
for label, jpath, dpath, rpath, zmin, zmax in CLASSES:
    J = json.load(open(jpath))
    prof = J["maps"]["lcdm"]["profile"]; fit = J["maps"]["lcdm"]["fit"]
    radii = np.array(prof["radii"], float)
    gam = np.array([np.nan if v is None else v for v in prof["gamma_star"]], float)
    n_used = len(load_cat(dpath, zmin, zmax)[2])  # same cut as the census run
    rra, rdec, rz, _ = load_cat(rpath, zmin, zmax)
    om = omega_sr(rra, rdec)
    d1, d2 = d_lcdm(np.array([zmin, zmax]))
    V = om/3.0*(d2**3 - d1**3)
    nbar = n_used/V
    m = np.isfinite(gam) & (gam > 0)
    rho = {p: float(np.interp(p, radii[m], gam[m])/nbar) for p in PIVOTS}
    rows.append({"class": label, "zwin": [zmin, zmax], "n_used": n_used,
                 "omega_sr": float(om), "V_Mpc3": float(V), "nbar": float(nbar),
                 "D2_lcdm": fit["D2"], "rho": rho})
    print(f"[{label}] z {zmin}-{zmax}  N={n_used:,}  Omega={om:.4f} sr  V={V:.3e} Mpc^3  nbar={nbar:.3e}")

lrg = next(r for r in rows if r["class"] == "LRG")
print("\nclass |   D2  | " + " | ".join(f"rho-1 @{int(p)}" for p in PIVOTS) + " | (bD)/(bD)_LRG @30,45")
for r in rows:
    ex = [r["rho"][p]-1.0 for p in PIVOTS]
    ratios = []
    for p in (30.0, 45.0):
        a, b = r["rho"][p]-1.0, lrg["rho"][p]-1.0
        ratios.append(np.sqrt(a/b) if a > 0 and b > 0 else float("nan"))
    print(f"{r['class']:5s} | {r['D2_lcdm']:.3f} | " + " | ".join(f"{v:9.4f}" for v in ex)
          + f" | {ratios[0]:.3f}, {ratios[1]:.3f}")
    r["bD_over_LRG_30"], r["bD_over_LRG_45"] = float(ratios[0]), float(ratios[1])

json.dump(rows, open("results/amplitude-check.json", "w"), indent=1)
print("\nwrote results/amplitude-check.json")
