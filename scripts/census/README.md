# P3 — The Census Fork Test: Runnable Pipeline

The conditional-census (P3) pipeline behind the Polyverse suite's census runs
(2026-08-23 → 2026-08-29), published per the corpus's pipeline-disclosure convention.
Companion page: [`the-polyverse-appendix.html`](../../the-polyverse-appendix.html) §A2
(run record); results narrative: [`the-polyverse-audit.html`](../../the-polyverse-audit.html) §9.2.

## What it computes

Conditional (Pietronero / Sylos Labini–Antal-style) integrated density
`Γ*(r) = ⟨n(<r)⟩`, averaged over **every embedded galaxy center**, on DESI DR1
LSS clustering catalogs, under **three z→distance maps**; the census dimension is
`D₂ = 3 + slope` of the weighted log–log fit.

- **Estimator** (`census_run.py::conditional_density`): per-center, per-scale
  embedding test against a **radially-resolved** random-count expectation
  (DESI randoms follow the survey n(z); a global mean density would be wrong by
  construction). Footprint solid angle from 0.5° angular cells of the randoms;
  radial profile n_ran(d) in 4 Mpc bins; keep center i at scale r iff
  `N_ran,i(<r) ≥ fill · E_i(r)` with `fill = 0.97` baseline (0.90 / 0.99 variants
  run as the embedding-criterion systematic). Weighted mean over kept centers
  (catalog `WEIGHT`), 100-fold bootstrap over centers, fixed seeds throughout.
- **Distance maps** (`DMAPS`): (a) flat ΛCDM H₀ = 70, Ωm = 0.3 — the
  Sylos Labini & Antal fiducial; (b) cosmography Taylor-3 and (c) Padé(2,1),
  both from **SN-only coefficients fixed in advance** (q₀ = −0.55, j₀ = 1.0,
  H₀ = 73; the overall scale cancels in slopes). The (b)–(c) spread is the
  truncation systematic. Firewall: no tuning; no FLRW assumed in (b)/(c).

## Data

- **DESI DR1 LSS clustering catalogs** (public, CC BY 4.0):
  https://data.desi.lbl.gov/doc/releases/dr1/ — tracers `BGS_BRIGHT-21.5`,
  `LRG`, `ELG_LOPnotqso`, `QSO`; both galactic caps; randoms
  `*_{0..}_clustering.ran.fits`. Files go under `data/` as
  `{TRACER}_{CAP}_clustering.dat.fits` + `{TRACER}_{CAP}_0_clustering.ran.fits`.
- **Mock catalogs** (public): AbacusSummit dark v4.2, root
  https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/
  — `complete` (fiber-complete), `ffa` (fast fiber assignment), and
  `altmtl0` (full alternative-MTL fiber assignment) flavors; `run_verdict*.sh`
  download the exact files.

## Run record (all runs on the gateway host, venv: numpy/scipy/astropy)

| date | script(s) | what ran |
|---|---|---|
| 2026-08-23 | `census_run.py` (v1 baseline) | 4 tracer classes × 2 caps, 3 distance maps; fits BGS 5–70, LRG/ELG 10–100, QSO 15–120 Mpc |
| 2026-08-27 | `run_hardening{,2,3}.sh` | z-slices within class; matched depth & scales; matched density (LRG→QSO thinning); 30–100 Mpc refit; fill 0.90/0.99; second random realization |
| 2026-08-28 | `run_hardening4.sh`, `run_hardening5.sh`, `qso_env_*.py`, `amplitude_check.py` | weight on/off (Δ = 0.000); catastrophic-z scramble f = 2/5/10/20%; conditional-amplitude cross-check; QSO environment test |
| 2026-08-29 | `run_verdict1.sh` … `run_verdict5.sh`, `mock_env_test.py`, `amplitude_mock.py` | wide-window turnover grids (to 260–300 Mpc, headline fits 100–260/300); AbacusSummit `complete` mock control through the identical pipeline (mock0 + realizations 1–3 for scatter); `ffa` leg; **full alt-MTL leg** (LRGalt0: D₂ = 2.922/2.906 vs data 2.921/2.907 — residual closed systematic); mock environment test |

Key outputs are kept small and included under `results/` (per-run JSON: full
Γ*(r) profiles, fits, settings echo; plus the four verdict-leg logs). The full
result set (66 JSONs) and downloaded catalogs (~GBs) live off-repo; every number
quoted on the suite pages traces to a JSON of this pipeline with settings inline.

## Reproduce

```bash
python3 -m venv venv && venv/bin/pip install numpy scipy astropy
# place catalogs under data/ (names above), then e.g.:
venv/bin/python census_run.py LRG NGC 0.4 1.1 10 100                 # v1 baseline
venv/bin/python census_run.py LRG NGC 0.4 1.1 100 260 --rgrid 10,260,26 --maps lcdm,pade21 --tag turnover260
./run_verdict2.sh   # mock-control leg (downloads AbacusSummit complete mock0)
./run_verdict5.sh   # alt-MTL leg (the LRG-residual decider)
```

`census_lcdm.py` is a memory-lean single-map wrapper. Flags:
`--ran N` (random realization), `--fill X` (embedding threshold), `--nmax N`
(matched-density subsample), `--scramble f` (catastrophic-z injection, n(z)
preserved exactly), `--noweight`, `--rgrid lo,hi,n`, `--maps a,b`, `--tag`.

Provenance: original working tree `~/.openclaw/workspace/p3-census/` (private);
published here 2026-09-07 as part of the Polyverse Phase-2 repairs. Nothing in
this folder is tuned to a target; coefficients and seeds are fixed in the source.
