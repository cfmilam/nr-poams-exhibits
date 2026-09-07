#!/usr/bin/env python3
"""
P3 — the census fork test.
Conditional (Pietronero/SL&A-style) integrated density Gamma*(r) = <n(<r)> from every
embedded galaxy center, on DESI DR1 clustering catalogs, under MULTIPLE z->distance maps:
  (a) fiducial flat LCDM  H0=70,  Om=0.3   [replicates Sylos Labini & Antal 2026 input]
  (b) cosmography Taylor-3 (SN-fit coefficients, no FLRW assumed)
  (c) cosmography Pade(2,1) (same coefficients; (b)-(c) spread = truncation systematic)
The question: does the census slope gamma [ <n(<r)> ~ r^-gamma, D2 = 3-gamma ] depend on
assuming FLRW geometry in the distance conversion? Nobody has published this. Firewall:
no tuning; coefficients fixed in advance from published SN-only cosmography.
"""
import sys, json, time
import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

C_KMS = 299792.458

# ---------------- distance maps ----------------
def d_lcdm(z, H0=70.0, Om=0.3, n=4096):
    """Comoving distance, flat LCDM (SL&A fiducial)."""
    zg = np.linspace(0, max(z.max(), 1.5), n)
    Ez = np.sqrt(Om*(1+zg)**3 + (1-Om))
    dc = np.concatenate([[0], np.cumsum(np.diff(zg)/((Ez[1:]+Ez[:-1])/2))]) * C_KMS/H0
    return np.interp(z, zg, dc)

# SN-only cosmography (fixed in advance; Pantheon+-era representative values):
# H0 = 73.0 (SN local calib; overall scale — cancels in slopes), q0 = -0.55, j0 = 1.0
Q0, J0, H0C = -0.55, 1.0, 73.0

def d_cosmo_taylor(z):
    """Luminosity-distance Taylor to 3rd order -> comoving d = d_L/(1+z). No FLRW."""
    dl = (C_KMS/H0C) * z * (1 + 0.5*(1-Q0)*z - (1/6)*(1 - Q0 - 3*Q0**2 + J0)*z**2)
    return dl/(1+z)

def d_cosmo_pade(z):
    """Pade(2,1) rearrangement of the same expansion (stable to z~1.1)."""
    a1 = 1.0
    a2 = 0.5*(1-Q0)
    a3 = -(1/6)*(1 - Q0 - 3*Q0**2 + J0)
    b1 = -a3/a2
    num = a1*z + (a2 + a1*b1)*z**2
    den = 1 + b1*z
    dl = (C_KMS/H0C) * num/den
    return dl/(1+z)

DMAPS = {"lcdm": d_lcdm, "taylor3": d_cosmo_taylor, "pade21": d_cosmo_pade}

# ---------------- IO ----------------
def load_cat(path, zmin, zmax):
    with fits.open(path) as h:
        d = h[1].data
        ra, dec, z = d["RA"].astype(np.float64), d["DEC"].astype(np.float64), d["Z"].astype(np.float64)
        w = (d["WEIGHT"] if "WEIGHT" in d.columns.names else np.ones_like(z)).astype(np.float64)
    m = (z > zmin) & (z < zmax) & np.isfinite(z)
    return ra[m], dec[m], z[m], w[m]

def to_cart(ra, dec, dist):
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    return np.column_stack([dist*np.cos(dec_r)*np.cos(ra_r),
                            dist*np.cos(dec_r)*np.sin(ra_r),
                            dist*np.sin(dec_r)])

# ---------------- estimator ----------------
def conditional_density(gal_xyz, gal_w, ran_xyz, radii, fill_thresh=0.97,
                        n_boot=100, rng=None):
    # fill_thresh is a v2-hardening knob (embedding-criterion systematic); default 0.97 = v1 baseline.
    """
    Gamma*(r) = weighted mean over embedded centers of N_i(<r) / (4/3 pi r^3).
    Embedding test per center per r: random-count fill fraction against a
    RADIALLY-RESOLVED expectation (DESI randoms follow the survey n(z); a global
    mean density is wrong by construction). For a sphere of radius r centred at
    radial distance d_i, the expected fully-embedded random count is
        E_i(r) = sum_d  n_ran(d) * pi (r^2 - (d - d_i)^2) * dd
    with n_ran(d) measured from the random catalogue itself in the SAME distance
    map (radial histogram / (Omega_sr * d^2 * dd)); Omega from coarse angular
    cells of the randoms. Self-contained: no fiducial cosmology re-enters.
    Keep centre i at scale r iff N_ran,i(<r) >= fill_thresh * E_i(r).
    Bootstrap over centers.
    """
    rng = rng or np.random.default_rng(42)
    tg, tr = cKDTree(gal_xyz), cKDTree(ran_xyz)
    # --- footprint solid angle from coarse angular cells of the randoms ---
    d_ran = np.linalg.norm(ran_xyz, axis=1)
    dec_r = np.arcsin(np.clip(ran_xyz[:, 2]/d_ran, -1, 1))
    ra_r = np.arctan2(ran_xyz[:, 1], ran_xyz[:, 0])
    cell = np.radians(0.5)                      # 0.5-deg cells
    iy = np.floor(dec_r/cell).astype(np.int64)
    ix = np.floor(ra_r*np.cos(dec_r)/cell).astype(np.int64)
    keys = ix*100000 + iy
    uniq, ucnt = np.unique(keys, return_counts=True)
    occ = uniq[ucnt >= 5]                       # occupied cells (noise floor)
    omega = len(occ)*cell*cell                  # steradians
    # --- radial density profile of randoms in this distance map ---
    dd = 4.0                                    # Mpc bins
    dmin, dmax = d_ran.min(), d_ran.max()
    edges = np.arange(dmin, dmax+dd, dd)
    hist, _ = np.histogram(d_ran, bins=edges)
    mid = 0.5*(edges[1:]+edges[:-1])
    n_prof = hist/(omega*mid**2*dd)             # n_ran(d) per Mpc^3
    n_of_d = lambda d: np.interp(d, mid, n_prof)
    # --- per-center radial distances ---
    d_gal = np.linalg.norm(gal_xyz, axis=1)
    out = {"radii": radii.tolist(), "gamma_star": [], "n_centers": [], "boot_err": []}
    for r in radii:
        vol = 4/3*np.pi*r**3
        # cheap radial pre-cut: sphere must fit inside the z-shell
        pre = (d_gal > dmin + r) & (d_gal < dmax - r)
        if pre.sum() < 50:
            out["gamma_star"].append(None); out["n_centers"].append(int(pre.sum())); out["boot_err"].append(None)
            continue
        # expected fully-embedded random count E_i(r), radially resolved
        slices = np.linspace(-r, r, 17)
        smid, sw = 0.5*(slices[1:]+slices[:-1]), np.diff(slices)
        E = np.zeros(pre.sum())
        dpre = d_gal[pre]
        for sm, w_ in zip(smid, sw):
            E += n_of_d(dpre + sm) * np.pi*(r*r - sm*sm) * w_
        cnt_ran = tr.query_ball_point(gal_xyz[pre], r, return_length=True)
        keep_sub = cnt_ran >= fill_thresh*E
        keep = np.zeros(len(gal_xyz), bool)
        keep[np.flatnonzero(pre)[keep_sub]] = True
        nk = int(keep.sum())
        if nk < 50:
            out["gamma_star"].append(None); out["n_centers"].append(nk); out["boot_err"].append(None)
            continue
        cnt = tg.query_ball_point(gal_xyz[keep], r, return_length=True).astype(np.float64) - 1.0  # exclude self
        dens = cnt/vol
        w = gal_w[keep]
        g = np.average(dens, weights=w)
        # bootstrap
        idx = np.arange(nk)
        bs = np.empty(n_boot)
        for b in range(n_boot):
            s = rng.choice(idx, nk, replace=True)
            bs[b] = np.average(dens[s], weights=w[s])
        out["gamma_star"].append(float(g))
        out["n_centers"].append(nk)
        out["boot_err"].append(float(bs.std()))
    return out

def fit_slope(radii, gam, err, rmin, rmax):
    r = np.asarray(radii, float); g = np.array([np.nan if v is None else v for v in gam], float)
    e = np.array([np.nan if v is None else v for v in err], float)
    m = (r >= rmin) & (r <= rmax) & np.isfinite(g) & np.isfinite(e) & (g > 0)
    if m.sum() < 3: return None
    x, y, s = np.log10(r[m]), np.log10(g[m]), e[m]/(g[m]*np.log(10))
    Wt = 1/s**2
    A = np.vstack([x, np.ones_like(x)]).T
    cov = np.linalg.inv(A.T @ (A*Wt[:, None]))
    coef = cov @ (A.T @ (y*Wt))
    slope, slope_err = coef[0], np.sqrt(cov[0, 0])
    return {"slope": float(slope), "slope_err": float(slope_err),
            "D2": float(3+slope), "rrange": [rmin, rmax], "npts": int(m.sum())}

# ---------------- main ----------------
def main(tracer, cap, zmin, zmax, rmin_fit, rmax_fit, ran_idx=0, fill=0.97, tag=None, nmax=None,
         scramble=None, noweight=False, rgrid=(2.0, 120.0, 22), maps=None):
    t0 = time.time()
    dat = f"data/{tracer}_{cap}_clustering.dat.fits"
    ran = f"data/{tracer}_{cap}_{ran_idx}_clustering.ran.fits"
    ra, dec, z, w = load_cat(dat, zmin, zmax)
    if nmax is not None and len(z) > nmax:
        # matched-density control: random subsample of galaxies (randoms untouched)
        sel = np.random.default_rng(7).choice(len(z), nmax, replace=False)
        ra, dec, z, w = ra[sel], dec[sel], z[sel], w[sel]
    if scramble:
        # catastrophic-z injection: permute the redshifts of a random fraction among
        # themselves (n(z) preserved EXACTLY; those galaxies' radial positions decorrelated
        # from their angular positions = in-window catastrophic-redshift model). Seed fixed.
        rng_s = np.random.default_rng(11)
        nsub = int(round(scramble*len(z)))
        idx = rng_s.choice(len(z), nsub, replace=False)
        z[idx] = z[rng_s.permutation(idx)]
        print(f"  [scramble] permuted z of {nsub:,}/{len(z):,} galaxies (f={scramble})", flush=True)
    if noweight:
        w = np.ones_like(w)
        print("  [noweight] WEIGHT column overridden to 1.0", flush=True)
    rra, rdec, rz, _ = load_cat(ran, zmin, zmax)
    print(f"[{tracer} {cap}] galaxies {len(z):,}  randoms {len(rz):,}  z in ({zmin},{zmax})  ran_idx={ran_idx} fill={fill} tag={tag}", flush=True)
    radii = np.geomspace(rgrid[0], rgrid[1], int(rgrid[2]))  # Mpc (comoving, in each map's units)
    results = {"tracer": tracer, "cap": cap, "zrange": [zmin, zmax], "ran_idx": ran_idx,
               "fill_thresh": fill, "tag": tag, "nmax": nmax, "scramble": scramble,
               "noweight": noweight, "n_gal_used": len(z), "rgrid": list(rgrid),
               "coeffs": {"Q0": Q0, "J0": J0, "H0_cosmo": H0C}, "maps": {}}
    for name, dmap in DMAPS.items():
        if maps and name not in maps:
            continue
        gx = to_cart(ra, dec, dmap(z))
        rx = to_cart(rra, rdec, dmap(rz))
        res = conditional_density(gx, w, rx, radii, fill_thresh=fill)
        fit = fit_slope(res["radii"], res["gamma_star"], res["boot_err"], rmin_fit, rmax_fit)
        results["maps"][name] = {"profile": res, "fit": fit}
        print(f"  [{name}] fit: {fit}", flush=True)
    suffix = f"_{tag}" if tag else ""
    out = f"results/{tracer}_{cap}{suffix}_census.json"
    json.dump(results, open(out, "w"), indent=1)
    print(f"wrote {out}  ({time.time()-t0:.0f}s)", flush=True)

if __name__ == "__main__":
    argv = sys.argv[1:]
    opts, skip, args = {}, set(), []
    for i, a in enumerate(argv):
        if a == "--ran":  opts["ran_idx"] = int(argv[i+1]);  skip.update({i, i+1})
        if a == "--fill": opts["fill"] = float(argv[i+1]);   skip.update({i, i+1})
        if a == "--tag":  opts["tag"] = argv[i+1];           skip.update({i, i+1})
        if a == "--nmax": opts["nmax"] = int(argv[i+1]);     skip.update({i, i+1})
        if a == "--scramble": opts["scramble"] = float(argv[i+1]); skip.update({i, i+1})
        if a == "--noweight": opts["noweight"] = True;        skip.add(i)
        if a == "--rgrid":
            lo, hi, n = argv[i+1].split(","); opts["rgrid"] = (float(lo), float(hi), int(n))
            skip.update({i, i+1})
        if a == "--maps": opts["maps"] = argv[i+1].split(","); skip.update({i, i+1})
    args = [a for i, a in enumerate(argv) if i not in skip]
    tracer, cap = args[0], args[1]
    zmin, zmax = float(args[2]), float(args[3])
    rmin, rmax = float(args[4]), float(args[5])
    main(tracer, cap, zmin, zmax, rmin, rmax, **opts)
