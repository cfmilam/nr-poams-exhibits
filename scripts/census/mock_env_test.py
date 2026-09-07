#!/usr/bin/env python3
"""
MOCK env test (2026-08-29) — verdict leg B3, closes the environment-interpretation question.
Same model-free measurement as qso_env_test.py stage (b), on AbacusSummit mock0 'complete':
(1+delta_env) around mock QSO vs mock LRG vs randoms, matched z 0.8-1.1, NGC, R=20 Mpc.
DATA result: QSO/LRG mean ratio 0.812, median 0.773 (~30sigma). If the LCDM mock (standard
HOD, no ledger-state selection) reproduces the under-dense QSO environments, the
ledger-state INTERPRETATION of the data split dies (standard halo occupancy suffices).
"""
import numpy as np
from scipy.spatial import cKDTree
from astropy.io import fits
from census_run import d_lcdm, to_cart

ZMIN, ZMAX = 0.8, 1.1
R_ENV = 20.0
RNG = np.random.default_rng(2029)

def load_xyz(path):
    with fits.open(path) as h:
        d = h[1].data
        ra, dec, z = d["RA"].astype(np.float64), d["DEC"].astype(np.float64), d["Z"].astype(np.float64)
    m = (z > ZMIN) & (z < ZMAX) & np.isfinite(z)
    return to_cart(ra[m], dec[m], d_lcdm(z[m]))

print("loading mock catalogs (NGC, z 0.8-1.1)...")
lrg = load_xyz("data/LRG_complete_NGC_clustering.dat.fits")
qso = load_xyz("data/QSO_complete_NGC_clustering.dat.fits")
ran = load_xyz("data/LRG_complete_NGC_0_clustering.ran.fits")
print(f"LRG {len(lrg):,}  QSO {len(qso):,}  RAN {len(ran):,}")

td, tr = cKDTree(lrg), cKDTree(ran)
Nd_tot, Nr_tot = float(len(lrg)), float(len(ran))

def delta_at(centers, self_in_data=False):
    Nd = td.query_ball_point(centers, R_ENV, return_length=True).astype(np.float64)
    Nr = tr.query_ball_point(centers, R_ENV, return_length=True).astype(np.float64)
    if self_in_data:
        Nd -= 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where(Nr > 0, (Nd/Nr)*(Nr_tot/Nd_tot) - 1.0, np.nan)
    return delta, Nr

def report(name, centers, self_in_data=False, nboot=200):
    delta, Nr = delta_at(centers, self_in_data)
    med = np.nanmedian(Nr)
    keep = np.isfinite(delta) & (Nr >= 0.5*med)   # edge trim, same as data run
    x = 1.0 + delta[keep]
    bs = np.empty(nboot)
    idx = np.arange(len(x))
    for b in range(nboot):
        bs[b] = x[RNG.choice(idx, len(x), replace=True)].mean()
    print(f"{name}: n={keep.sum():,}  mean={x.mean():.3f}+/-{bs.std():.3f}  "
          f"median={np.median(x):.3f}  p10={np.percentile(x,10):.2f} p90={np.percentile(x,90):.2f}")
    return x.mean(), np.median(x)

# randoms self-sample (sanity ~1.0): subsample for speed
sub = ran[RNG.choice(len(ran), 300_000, replace=False)]
m_r, md_r = report("RAND", sub)
m_l, md_l = report("LRG ", lrg, self_in_data=True)
m_q, md_q = report("QSO ", qso)

print(f"\nMOCK QSO/LRG ratio: mean {m_q/m_l:.3f}  median {md_q/md_l:.3f}")
print("DATA QSO/LRG ratio: mean 0.812  median 0.773  (R=20 Mpc)")
