#!/usr/bin/env python3
"""Amplitude + D2 extraction for the reweighted-LRG gamma grid.
conditional amplitude metric (matches amplitude_check.py): (bD)/(bD)_LRG = sqrt((rho_x-1)/(rho_LRG-1))
rho = Gamma*/nbar. For weighted runs nbar = sum(WEIGHT in-window)/V (same Omega,V as LRG baseline).
Target from data: QSO D2=2.795, QSO (bD)/(bD)_LRG = 0.55 (scale-flat 30-70)."""
import json, glob, numpy as np
from astropy.io import fits
from census_run import d_lcdm, load_cat
from amplitude_check import omega_sr

ZMIN, ZMAX = 0.8, 1.1
PIV = [30.0, 45.0, 70.0]

def nbar_from(dat, ran):
    with fits.open(dat) as h:
        d = h[1].data; z = d["Z"].astype(float); w = d["WEIGHT"].astype(float)
    m = (z > ZMIN) & (z < ZMAX) & np.isfinite(z)
    sw = w[m].sum()
    rra, rdec, rz, _ = load_cat(ran, ZMIN, ZMAX)
    om = omega_sr(rra, rdec)
    d1, d2 = d_lcdm(np.array([ZMIN, ZMAX]))
    V = om / 3.0 * (d2**3 - d1**3)
    return sw / V, sw, V

def rho_minus1(jpath, nbar):
    J = json.load(open(jpath))
    prof = J["maps"]["lcdm"]["profile"]; fit = J["maps"]["lcdm"]["fit"]
    r = np.array(prof["radii"], float)
    g = np.array([np.nan if v is None else v for v in prof["gamma_star"]], float)
    m = np.isfinite(g) & (g > 0)
    rho = {p: float(np.interp(p, r[m], g[m]) / nbar) for p in PIV}
    return fit["D2"], {p: rho[p] - 1.0 for p in PIV}

# LRG unweighted baseline
lrg_nbar, _, _ = nbar_from("data/LRG_NGC_clustering.dat.fits", "data/LRG_NGC_0_clustering.ran.fits")
D2_lrg, rm1_lrg = rho_minus1("results/LRG_NGC_z08-11_census.json", lrg_nbar)

rows = []
GAMS = [("g0p5", 0.5), ("g1p0", 1.0), ("g1p5", 1.5), ("g2p0", 2.0), ("g3p0", 3.0)]
for tagg, g in GAMS:
    dat = f"data/LRGenvR20{tagg}_NGC_clustering.dat.fits"
    ran = f"data/LRGenvR20{tagg}_NGC_0_clustering.ran.fits"
    j = f"results/LRGenvR20{tagg}_NGC_env_census.json"
    nb, sw, V = nbar_from(dat, ran)
    D2, rm1 = rho_minus1(j, nb)
    bD = {p: float(np.sqrt(rm1[p] / rm1_lrg[p])) if rm1[p] > 0 and rm1_lrg[p] > 0 else float("nan") for p in PIV}
    bD_mean = float(np.nanmean([bD[p] for p in PIV]))
    rows.append({"gamma": g, "D2": D2, "bD_over_LRG": bD, "bD_mean": bD_mean})
    print(f"gamma={g:<4} D2={D2:.3f}  (bD)/(bD)_LRG @30/45/70 = {bD[30.0]:.3f}/{bD[45.0]:.3f}/{bD[70.0]:.3f} (mean {bD_mean:.3f})")

print(f"\nLRG baseline: D2={D2_lrg:.3f}  rho-1@30/45/70={rm1_lrg[30.0]:.3f}/{rm1_lrg[45.0]:.3f}/{rm1_lrg[70.0]:.3f}")
print("TARGET (QSO): D2=2.795  (bD)/(bD)_LRG=0.55")

# interpolate gamma that hits each target
gs = np.array([r["gamma"] for r in rows])
d2s = np.array([r["D2"] for r in rows])
amps = np.array([r["bD_mean"] for r in rows])
def solve(xs, ys, target):
    o = np.argsort(ys)
    return float(np.interp(target, ys[o], xs[o]))
g_D2 = solve(gs, d2s, 2.795)
g_amp = solve(gs, amps, 0.55)
print(f"\ngamma_D2  (hits D2=2.795)  = {g_D2:.2f}")
print(f"gamma_amp (hits amp=0.55) = {g_amp:.2f}")
json.dump({"rows": rows, "D2_lrg": D2_lrg, "rm1_lrg": rm1_lrg,
           "gamma_D2": g_D2, "gamma_amp": g_amp,
           "target": {"D2": 2.795, "amp": 0.55}},
          open("results/qso-env-model.json", "w"), indent=1)
print("wrote results/qso-env-model.json")
