#!/bin/bash
# P3 VERDICT LEG B — LCDM mock control through the IDENTICAL pipeline (2026-08-29)
# AbacusSummit DR1 mock0 (Planck-2018 flat LCDM, Om~0.3137), 'complete' flavor (no FA loss).
# Questions: (1) does a LCDM universe processed by OUR estimator homogenize (D2(r)->3
# at large r)? (2) does mock QSO sit at its bias station (data QSO reads +0.27 above)?
# Caveat noted: z->d maps identical to data runs (SL&A fiducial H0=70 Om=0.30 vs mock
# truth Om=0.3137) — map-spread tests bound this at <=0.05 in D2; identical-pipe is the point.
set -e
cd "$(dirname "$0")"
PY=venv/bin/python
BASE=https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/AbacusSummit/dark/v4.2/mock0
LOG=results/verdict-mock-20260829.log
{
echo "=== VERDICT LEG B launched $(date) ==="
for f in LRG_complete_NGC_clustering.dat.fits LRG_complete_NGC_0_clustering.ran.fits \
         QSO_complete_NGC_clustering.dat.fits QSO_complete_NGC_0_clustering.ran.fits; do
  if [ ! -s "data/$f" ]; then
    echo "--- downloading $f"
    curl -fSL --retry 5 --retry-delay 10 -C - -o "data/$f" "$BASE/$f"
  else
    echo "--- have $f"
  fi
done
ls -l data/*_complete_* || true
echo '--- MOCK LRG full window 0.4-1.1, baseline fit 10-100 (data: 2.517) ---'
$PY census_run.py LRG_complete NGC 0.4 1.1 10 100 --maps lcdm --tag mock0
echo '--- MOCK LRG matched-z 0.8-1.1, fit 15-100 (data: 2.529) ---'
$PY census_run.py LRG_complete NGC 0.8 1.1 15 100 --maps lcdm --tag mock0-z08-11
echo '--- MOCK QSO matched-z 0.8-1.1, fit 15-100 (data: 2.795 = the anomaly) ---'
$PY census_run.py QSO_complete NGC 0.8 1.1 15 100 --maps lcdm --tag mock0-z08-11
echo '--- MOCK QSO full 0.8-1.8, fit 15-120 (data: 2.828) ---'
$PY census_run.py QSO_complete NGC 0.8 1.8 15 120 --maps lcdm --tag mock0-full
echo '--- MOCK LRG TURNOVER grid 10-260, fit 100-260 (THE control: does LCDM homogenize in our pipe?) ---'
$PY census_run.py LRG_complete NGC 0.4 1.1 100 260 --rgrid 10,260,26 --maps lcdm --tag mock0-turnover260
echo '--- MOCK QSO TURNOVER grid 15-300, fit 120-300 ---'
$PY census_run.py QSO_complete NGC 0.8 1.8 120 300 --rgrid 15,300,24 --maps lcdm --tag mock0-turnover300
echo "=== VERDICT LEG B complete $(date) ==="
} >> "$LOG" 2>&1
