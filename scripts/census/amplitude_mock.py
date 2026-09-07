#!/usr/bin/env python3
"""
MOCK amplitude check (2026-08-29) — verdict leg B2.
Same internal amplitude statistic as amplitude_check.py, on AbacusSummit mock0
'complete' LRG + QSO at matched z 0.8-1.1 (lcdm map, NGC), vs the DATA values.
Question: does the LCDM mock reproduce QSO's conditional-amplitude suppression
(data: QSO rho-1 = 0.55x LRG, scale-flat 30-70), or does the amplitude leg of the
QSO anomaly survive the mock control?
"""
import json, numpy as np
from astropy.io import fits
from census_run import d_lcdm, load_cat

C = 0.5
def omega_sr(ra, dec):
    cell = np.radians(C)
    dec_r, ra_r = np.radians(dec), np.radians(ra)
    iy = np.floor(dec_r/cell).astype(np.int64)
    ix = np.floor(ra_r*np.cos(dec_r)/cell).astype(np.int64)
    keys = ix*100000 + iy
    uniq, ucnt = np.unique(keys, return_counts=True)
    return (ucnt >= 5).sum()*cell*cell

CLASSES = [
    ("mockLRG", "results/LRG_complete_NGC_mock0-z08-11_census.json",
     "data/LRG_complete_NGC_clustering.dat.fits", "data/LRG_complete_NGC_0_clustering.ran.fits", 0.8, 1.1),
    ("mockQSO", "results/QSO_complete_NGC_mock0-z08-11_census.json",
     "data/QSO_complete_NGC_clustering.dat.fits", "data/QSO_complete_NGC_0_clustering.ran.fits", 0.8, 1.1),
]
PIVOTS = [20.0, 30.0, 45.0, 70.0]

out = {}
for label, jpath, dpath, rpath, zmin, zmax in CLASSES:
    r = json.load(open(jpath))
    prof = r["maps"]["lcdm"]["profile"]
    radii = np.array(prof["radii"], float)
    gam = np.array([x if x is not None else np.nan for x in prof["gamma_star"]], float)
    ra, dec, z, w = load_cat(dpath, zmin, zmax)
    rra, rdec, rz, _ = load_cat(rpath, zmin, zmax)
    omega = omega_sr(rra, rdec)
    d1, d2 = d_lcdm(np.array([zmin])), d_lcdm(np.array([zmax]))
    V = omega/3.0*(d2[0]**3 - d1[0]**3)
    nbar = len(z)/V
    rho = gam/nbar
    out[label] = {}
    print(f"{label}: N={len(z):,} Omega={omega:.3f} sr V={V:.3e} nbar={nbar:.3e}")
    for p in PIVOTS:
        i = int(np.argmin(np.abs(radii-p)))
        out[label][p] = rho[i]-1.0
        print(f"   rho-1 @ {radii[i]:.0f} Mpc = {rho[i]-1:.3f}")

print("\n=== MOCK QSO/LRG conditional-amplitude ratio (data: 0.545/0.547/0.548 @30/45/70) ===")
for p in PIVOTS:
    r_ = out["mockQSO"][p]/out["mockLRG"][p]
    print(f"   ratio @ {p:.0f} Mpc = {r_:.3f}")
