#!/bin/bash
# ROUND 2 leg 2 — FULL alt-MTL fiber-assignment mock (most data-like) for the LRG residual (2026-08-29)
# ffa closed ~29%/~48% of the data-mock gap; altmtl decides how much of the rest is survey plumbing.
set -e
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/verdict-altmtl-20260829.log
{
echo "=== ROUND2 ALTMTL LEG launched $(date) ==="
BASE=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/AbacusSummit/dark/v4.2/altmtl0/mock0/LSScats
echo "--- discovering catalog names at $BASE/"
IDX=$(curl -fsSL "$BASE/" || true)
DAT=$(echo "$IDX" | grep -o 'href="[^"]*LRG[^"]*NGC[^"]*clustering\.dat\.fits"' | head -1 | sed 's/href="//;s/"//')
RAN=$(echo "$IDX" | grep -o 'href="[^"]*LRG[^"]*NGC[^"]*0_clustering\.ran\.fits"' | head -1 | sed 's/href="//;s/"//')
echo "found DAT=$DAT RAN=$RAN"
if [ -z "$DAT" ] || [ -z "$RAN" ]; then echo "DISCOVERY FAILED — index listing follows:"; echo "$IDX" | grep -o 'href="[^"]*"' | head -40; exit 1; fi
for f in "$DAT" "$RAN"; do
  OUT="data/altmtl0-$f"
  if [ ! -s "$OUT" ]; then
    echo "--- downloading altmtl0 $f"
    curl -fSL --retry 5 --retry-delay 10 -C - -o "$OUT" "$BASE/$f"
  fi
done
ln -sf "altmtl0-$DAT" "data/LRGalt0_NGC_clustering.dat.fits"
ln -sf "altmtl0-$RAN" "data/LRGalt0_NGC_0_clustering.ran.fits"
echo "--- ALTMTL0 LRG turnover 0.4-1.1, grid 10-260, fit 100-260 ---"
$PY census_run.py LRGalt0 NGC 0.4 1.1 100 260 --rgrid 10,260,26 --maps lcdm --tag turnover260
echo "=== ROUND2 ALTMTL LEG complete $(date) ==="
$PY - << 'PYEOF'
import json, numpy as np
r = json.load(open("results/LRGalt0_NGC_turnover260_census.json"))
p = r["maps"]["lcdm"]["profile"]
rr = np.array(p["radii"], float)
g = np.array([x if x is not None else np.nan for x in p["gamma_star"]], float)
ok = np.isfinite(g) & (g > 0)
lr, lg = np.log(rr[ok]), np.log(g[ok])
print("ALTMTL0 windows (complete 2.962/2.965 | ffa 2.950/2.937 | DATA 2.921/2.907):")
for lo, hi in [(150, 220), (200, 300)]:
    m2 = (rr[ok] >= lo) & (rr[ok] <= hi)
    d2 = 3 + np.polyfit(lr[m2], lg[m2], 1)[0]
    print(f"  {lo}-{hi}: D2 = {d2:.3f}")
print("Guide: altmtl ≈ 2.92 → survey plumbing owns the residual (anomaly CLOSED).")
print("       altmtl ≈ 2.94+ → remainder survives full fiber treatment; imaging weights next.")
PYEOF
} >> "$LOG" 2>&1
