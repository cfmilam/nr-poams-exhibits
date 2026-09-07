#!/usr/bin/env python3
"""
P3 PHYSICAL-READING TEST — QSO ledger-state selection, one-parameter over-constrained.

Hypothesis: QSO activity selects seats whose local books are still being written; such
seats preferentially sit in LESS-settled (lower-overdensity) large-scale environments.
Model:  P(active|env) prop (1+delta_env)^(-gamma), ONE parameter gamma.
It must reproduce BOTH standing QSO observables from one number:
   (i)  D2 = 2.795  (vs LRG 2.529; gap +0.27)
   (ii) conditional amplitude (bD)/(bD)_LRG = 0.55  (scale-flat 30-70)
The gamma that reproduces (i) must also reproduce (ii), or the model dies.

Density tracer = LRG catalog (matched-z 0.8-1.1). Environment field via KDTree top-hat.
Stages:
  (a) build delta_env at R_env for QSO, LRG, and random-point positions
  (b) MEASURE (model-free): (1+delta_env) distributions QSO vs LRG vs random
  (c) write reweighted LRG dat.fits w_i = W_i*(1+delta_i)^(-gamma) for a gamma grid
      (census_run.py then re-runs the conditional census off the WEIGHT column)
"""
import sys, json, time
import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
from census_run import d_lcdm, to_cart, load_cat

ZMIN, ZMAX = 0.8, 1.1
RNG = np.random.default_rng(2029)

def load_xyz(path, zmin=ZMIN, zmax=ZMAX, cols_w=True):
    with fits.open(path) as h:
        d = h[1].data
        ra, dec, z = d["RA"].astype(np.float64), d["DEC"].astype(np.float64), d["Z"].astype(np.float64)
        w = (d["WEIGHT"] if ("WEIGHT" in d.columns.names and cols_w) else np.ones_like(z)).astype(np.float64)
    m = (z > zmin) & (z < zmax) & np.isfinite(z)
    ra, dec, z, w = ra[m], dec[m], z[m], w[m]
    xyz = to_cart(ra, dec, d_lcdm(z))
    return xyz, w, m

def build_field(cap):
    """LRG data + randoms (this cap) -> KDTrees for the density field."""
    lrg_xyz, lrg_w, _ = load_xyz(f"data/LRG_{cap}_clustering.dat.fits")
    ran_xyz, _, _ = load_xyz(f"data/LRG_{cap}_0_clustering.ran.fits")
    return lrg_xyz, lrg_w, ran_xyz

def delta_at(centers, td, tr, Nd_tot, Nr_tot, R, self_in_data=False):
    """delta_env = (Nd/Nr)*(Nr_tot/Nd_tot) - 1 in top-hat sphere R (randoms handle
    survey window + n(z) selection + edges). self_in_data: subtract 1 from Nd (center is
    itself an LRG). Returns delta and the raw random count Nr (for edge trimming)."""
    Nd = td.query_ball_point(centers, R, return_length=True).astype(np.float64)
    Nr = tr.query_ball_point(centers, R, return_length=True).astype(np.float64)
    if self_in_data:
        Nd = Nd - 1.0
    scale = Nr_tot / Nd_tot
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where(Nr > 0, (Nd / Nr) * scale - 1.0, np.nan)
    return delta, Nr

def pct_block(name, x, nboot=200):
    x = x[np.isfinite(x)]
    ps = [5, 10, 25, 50, 75, 90, 95]
    perc = np.percentile(x, ps)
    mean = x.mean()
    # bootstrap error on the mean and median
    n = len(x)
    bm = np.empty(nboot); bmed = np.empty(nboot)
    for b in range(nboot):
        s = x[RNG.integers(0, n, n)]
        bm[b] = s.mean(); bmed[b] = np.median(s)
    return {"name": name, "n": int(n), "mean": float(mean), "mean_err": float(bm.std()),
            "median": float(perc[3]), "median_err": float(bmed.std()),
            "pcts": {str(p): float(v) for p, v in zip(ps, perc)}}

def measure(R, log):
    """Stage (b): model-free (1+delta) distributions QSO vs LRG vs random, both caps pooled."""
    onep = {"QSO": [], "LRG": [], "RAND": []}
    for cap in ("NGC", "SGC"):
        lrg_xyz, lrg_w, ran_xyz = build_field(cap)
        qso_xyz, _, _ = load_xyz(f"data/QSO_{cap}_clustering.dat.fits")
        td, tr = cKDTree(lrg_xyz), cKDTree(ran_xyz)
        Nd_tot, Nr_tot = len(lrg_xyz), len(ran_xyz)
        # random-point sample = subsample of the randoms themselves
        rsub = ran_xyz[RNG.choice(Nr_tot, min(200_000, Nr_tot), replace=False)]
        # LRG centers subsample for speed (distribution only)
        lsub_idx = RNG.choice(Nd_tot, min(200_000, Nd_tot), replace=False)
        dq, nrq = delta_at(qso_xyz, td, tr, Nd_tot, Nr_tot, R)
        dl, nrl = delta_at(lrg_xyz[lsub_idx], td, tr, Nd_tot, Nr_tot, R, self_in_data=True)
        dr, nrr = delta_at(rsub, td, tr, Nd_tot, Nr_tot, R)
        # edge trim: common threshold from the random sample's own Nr (interior tracer)
        thr = 0.5 * np.median(nrr[np.isfinite(dr)])
        onep["QSO"].append((1 + dq)[(nrq >= thr) & np.isfinite(dq)])
        onep["LRG"].append((1 + dl)[(nrl >= thr) & np.isfinite(dl)])
        onep["RAND"].append((1 + dr)[(nrr >= thr) & np.isfinite(dr)])
        print(f"  [{cap}] R={R}  Nd={Nd_tot:,} Nr={Nr_tot:,}  edge thr(Nr)>={thr:.0f}  "
              f"QSO {len(onep['QSO'][-1]):,} LRG {len(onep['LRG'][-1]):,} RAND {len(onep['RAND'][-1]):,}", flush=True)
        log.write(f"[measure {cap} R={R}] thr={thr:.0f} keptQ={len(onep['QSO'][-1])} keptL={len(onep['LRG'][-1])} keptR={len(onep['RAND'][-1])}\n")
    res = {k: pct_block(k, np.concatenate(v)) for k, v in onep.items()}
    return res

def write_reweighted(cap, R, gammas, log):
    """Stage (c): delta_env for ALL LRG data in-window (this cap), write reweighted dat.fits
    per gamma. census_run reads WEIGHT; centers with bad delta (edge) get weight 0 -> excluded
    from the weighted census average (embedding still uses randoms, unaffected)."""
    lrg_xyz, lrg_w, ran_xyz = build_field(cap)
    td, tr = cKDTree(lrg_xyz), cKDTree(ran_xyz)
    Nd_tot, Nr_tot = len(lrg_xyz), len(ran_xyz)
    delta, Nr = delta_at(lrg_xyz, td, tr, Nd_tot, Nr_tot, R, self_in_data=True)
    thr = 0.5 * np.median(Nr[np.isfinite(delta)])
    good = np.isfinite(delta) & (Nr >= thr) & (delta > -1.0)
    onep = np.clip(1.0 + delta, 1e-6, None)
    # need the in-window mask to align rows with the source fits WEIGHT column
    src = f"data/LRG_{cap}_clustering.dat.fits"
    with fits.open(src) as h:
        d = h[1].data
        z = d["Z"].astype(np.float64)
        inwin = (z > ZMIN) & (z < ZMAX) & np.isfinite(z)
    print(f"  [reweight {cap}] R={R} LRG in-window {inwin.sum():,}  good delta {good.sum():,}  "
          f"median(1+delta)={np.median(onep[good]):.3f}", flush=True)
    log.write(f"[reweight {cap} R={R}] inwin={int(inwin.sum())} good={int(good.sum())} "
              f"med(1+d)={np.median(onep[good]):.4f} thr={thr:.0f}\n")
    paths = {}
    for g in gammas:
        w_new = np.where(good, lrg_w * onep ** (-g), 0.0)
        # write a full-length WEIGHT array back to a copy of the fits (only WEIGHT changed)
        with fits.open(src) as h:
            cols = h[1].columns
            data = h[1].data
            wfull = np.array(data["WEIGHT"], dtype=np.float64)
            # map in-window new weights into full array; out-of-window rows keep orig (census cuts them)
            wfull_inwin = wfull.copy()
            wfull_inwin[inwin] = w_new
            newcols = []
            for c in cols:
                if c.name == "WEIGHT":
                    newcols.append(fits.Column(name="WEIGHT", format=c.format, array=wfull_inwin))
                else:
                    newcols.append(fits.Column(name=c.name, format=c.format, array=data[c.name]))
            hdu = fits.BinTableHDU.from_columns(newcols)
        tag = f"envR{int(R)}g{str(g).replace('.','p')}"
        outdat = f"data/LRG_{cap}_{tag}_clustering.dat.fits"
        hdu.writeto(outdat, overwrite=True)
        # census_run expects data/{tracer}_{cap}_clustering.dat.fits with tracer=LRG_..{tag}?
        # simpler: we register a synthetic tracer name by symlinking ran to match.
        import os
        ransrc = f"LRG_{cap}_0_clustering.ran.fits"
        ranlink = f"data/LRG{tag}_{cap}_0_clustering.ran.fits"
        datlink = f"data/LRG{tag}_{cap}_clustering.dat.fits"
        for lk, tgt in ((ranlink, ransrc), (datlink, f"LRG_{cap}_{tag}_clustering.dat.fits")):
            if os.path.islink(lk) or os.path.exists(lk):
                os.remove(lk)
            os.symlink(tgt, lk)
        paths[g] = (f"LRG{tag}", outdat)
    return paths, {"good": int(good.sum()), "inwin": int(inwin.sum()),
                   "median_1pd": float(np.median(onep[good]))}

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "measure"
    logpath = "results/qso-env-20260829.log"
    log = open(logpath, "a")
    log.write(f"\n===== qso_env_test {mode} {time.strftime('%Y-%m-%d %H:%M')} =====\n")
    if mode == "measure":
        out = {}
        for R in (20.0, 40.0):
            print(f"[MEASURE R={R}]", flush=True)
            out[str(R)] = measure(R, log)
        json.dump(out, open("results/qso-env-measure.json", "w"), indent=1)
        for R, res in out.items():
            log.write(f"\n-- R={R} --\n")
            for k in ("RAND", "LRG", "QSO"):
                b = res[k]
                log.write(f"{k}: mean {b['mean']:.4f}+/-{b['mean_err']:.4f}  median {b['median']:.4f}+/-{b['median_err']:.4f}  "
                          f"p10={b['pcts']['10']:.3f} p50={b['pcts']['50']:.3f} p90={b['pcts']['90']:.3f}\n")
        print("wrote results/qso-env-measure.json")
    elif mode == "reweight":
        R = float(sys.argv[2]); cap = sys.argv[3] if len(sys.argv) > 3 else "NGC"
        gammas = [0.5, 1.0, 1.5, 2.0, 3.0]
        paths, meta = write_reweighted(cap, R, gammas, log)
        json.dump({"paths": {str(g): p for g, p in paths.items()}, "meta": meta},
                  open(f"results/qso-env-reweight-{cap}-R{int(R)}.json", "w"), indent=1)
        print("tracers:", {g: p[0] for g, p in paths.items()})
    log.close()
