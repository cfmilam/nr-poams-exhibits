#!/bin/bash
# P3 VERDICT LEG C — mock realization scatter for the LRG large-scale residual (2026-08-29)
# Data LRG plateaus at 2.91 (150-300 Mpc); mock0 climbs to 2.96-2.97. Is Δ≈0.05 real?
# Kill test: run mock1-3 LRG 'complete' turnovers; realization scatter decides.
set -e
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/verdict-mockscatter-20260829.log
{
echo "=== VERDICT LEG C launched $(date) ==="
for i in 1 2 3; do
  BASE=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/AbacusSummit/dark/v4.2/mock$i
  for f in LRG_complete_NGC_clustering.dat.fits LRG_complete_NGC_0_clustering.ran.fits; do
    OUT="data/mock$i-$f"
    if [ ! -s "$OUT" ]; then
      echo "--- downloading mock$i $f"
      curl -fSL --retry 5 --retry-delay 10 -C - -o "$OUT" "$BASE/$f"
    fi
  done
  # census_run.py expects data/{tracer}_{cap}_clustering.dat.fits naming; symlink per-mock
  ln -sf "mock${i}-LRG_complete_NGC_clustering.dat.fits" "data/LRGmock${i}_NGC_clustering.dat.fits"
  ln -sf "mock${i}-LRG_complete_NGC_0_clustering.ran.fits" "data/LRGmock${i}_NGC_0_clustering.ran.fits"
  echo "--- MOCK$i LRG turnover 0.4-1.1, grid 10-260, fit 100-260 ---"
  $PY census_run.py LRGmock$i NGC 0.4 1.1 100 260 --rgrid 10,260,26 --maps lcdm --tag turnover260
done
echo "=== VERDICT LEG C complete $(date) ==="
echo "Scatter summary:"
$PY - << 'PYEOF'
import json, numpy as np
vals = {}
for tag in ["LRG_complete_NGC_mock0-turnover260", "LRGmock1_NGC_turnover260",
            "LRGmock2_NGC_turnover260", "LRGmock3_NGC_turnover260"]:
    try:
        r = json.load(open(f"results/{tag}_census.json"))
        p = r["maps"]["lcdm"]["profile"]
        rr = np.array(p["radii"], float)
        g = np.array([x if x is not None else np.nan for x in p["gamma_star"]], float)
        ok = np.isfinite(g) & (g > 0)
        lr, lg = np.log(rr[ok]), np.log(g[ok])
        out = []
        for lo, hi in [(150, 220), (200, 300)]:
            m2 = (rr[ok] >= lo) & (rr[ok] <= hi)
            out.append(3 + np.polyfit(lr[m2], lg[m2], 1)[0] if m2.sum() >= 3 else np.nan)
        vals[tag] = out
        print(tag, [f"{v:.3f}" for v in out])
    except Exception as e:
        print(tag, "missing", e)
arr = np.array([v for v in vals.values() if np.all(np.isfinite(v))])
if len(arr) >= 3:
    print(f"mock mean 150-220: {arr[:,0].mean():.3f} +/- {arr[:,0].std():.3f}   (DATA 2.921)")
    print(f"mock mean 200-300: {arr[:,1].mean():.3f} +/- {arr[:,1].std():.3f}   (DATA 2.907)")
    print(f"pull 150-220: {(arr[:,0].mean()-2.921)/max(arr[:,0].std(),1e-9):.1f} sigma_realization")
    print(f"pull 200-300: {(arr[:,1].mean()-2.907)/max(arr[:,1].std(),1e-9):.1f} sigma_realization")
PYEOF
} >> "$LOG" 2>&1
