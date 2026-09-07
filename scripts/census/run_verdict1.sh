#!/bin/bash
# P3 VERDICT LEG A — homogenization turnover search (2026-08-29)
# Question: does the running conditional dimension D2(r) approach 3 at the largest
# embeddable scales (the homogenization verdict), or stay at its class value?
# Full Gamma*(r) profiles to 150-300 Mpc; headline fit = large-scale window only.
# Maps: lcdm + pade21 (taylor3 known invalid past z~1.5; pade is the model-independent leg).
set -e
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/verdict-turnover-20260829.log
{
echo "=== VERDICT LEG A launched $(date) ==="
echo '--- LRG NGC 0.4-1.1, grid 10-260, headline fit 100-260 ---'
$PY census_run.py LRG NGC 0.4 1.1 100 260 --rgrid 10,260,26 --maps lcdm,pade21 --tag turnover260
echo '--- ELG NGC 0.8-1.6, grid 10-300, headline fit 100-300 ---'
$PY census_run.py ELG_LOPnotqso NGC 0.8 1.6 100 300 --rgrid 10,300,26 --maps lcdm,pade21 --tag turnover300
echo '--- BGS NGC 0.05-0.40, grid 5-150, headline fit 70-150 ---'
$PY census_run.py BGS_BRIGHT-21.5 NGC 0.05 0.40 70 150 --rgrid 5,150,24 --maps lcdm,pade21 --tag turnover150
echo '--- QSO NGC 0.8-1.8, grid 15-300, headline fit 120-300 ---'
$PY census_run.py QSO NGC 0.8 1.8 120 300 --rgrid 15,300,24 --maps lcdm,pade21 --tag turnover300
echo "=== VERDICT LEG A complete $(date) ==="
} >> "$LOG" 2>&1
