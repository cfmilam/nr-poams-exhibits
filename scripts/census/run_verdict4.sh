#!/bin/bash
# ROUND 2 (anomaly systematics leg) — fiber-assignment kill test for the LRG plateau residual (2026-08-29)
# Data LRG plateaus 2.907-2.921 (150-300 Mpc); complete mocks 2.962-2.965 (±0.002 realization).
# Suspect: data carries fiber-assignment incompleteness; complete mocks don't.
# Test: identical pipeline on AbacusSummit mock0 LRG *ffa* (fast-fiber-assign) NGC.
# If ffa drops D2 by ~0.04 → systematics owns the residual. If it stays ~2.96 → anomaly survives.
set -e
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/verdict-ffa-20260829.log
{
echo "=== ROUND2 FFA LEG launched $(date) ==="
BASE=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/AbacusSummit/dark/v4.2/mock0
for f in LRG_ffa_NGC_clustering.dat.fits LRG_ffa_NGC_0_clustering.ran.fits; do
  OUT="data/mock0-$f"
  if [ ! -s "$OUT" ]; then
    echo "--- downloading mock0 $f"
    curl -fSL --retry 5 --retry-delay 10 -C - -o "$OUT" "$BASE/$f"
  fi
done
ln -sf "mock0-LRG_ffa_NGC_clustering.dat.fits" "data/LRGffa0_NGC_clustering.dat.fits"
ln -sf "mock0-LRG_ffa_NGC_0_clustering.ran.fits" "data/LRGffa0_NGC_0_clustering.ran.fits"
echo "--- MOCK0 LRG FFA turnover 0.4-1.1, grid 10-260, fit 100-260 ---"
$PY census_run.py LRGffa0 NGC 0.4 1.1 100 260 --rgrid 10,260,26 --maps lcdm --tag turnover260
echo "=== ROUND2 FFA LEG complete $(date) ==="
$PY - << 'PYEOF'
import json, numpy as np
r = json.load(open("results/LRGffa0_NGC_turnover260_census.json"))
p = r["maps"]["lcdm"]["profile"]
rr = np.array(p["radii"], float)
g = np.array([x if x is not None else np.nan for x in p["gamma_star"]], float)
ok = np.isfinite(g) & (g > 0)
lr, lg = np.log(rr[ok]), np.log(g[ok])
print("FFA mock0 windows:")
for lo, hi in [(150, 220), (200, 300)]:
    m2 = (rr[ok] >= lo) & (rr[ok] <= hi)
    d2 = 3 + np.polyfit(lr[m2], lg[m2], 1)[0]
    print(f"  {lo}-{hi}: D2 = {d2:.3f}   (complete mocks 2.962/2.965 ± 0.002; DATA 2.921/2.907)")
print("Verdict guide: ffa ≈ 2.92 → fiber assignment owns the residual (anomaly dies).")
print("               ffa ≈ 2.96 → residual survives fiber assignment; next: imaging weights.")
PYEOF
} >> "$LOG" 2>&1
