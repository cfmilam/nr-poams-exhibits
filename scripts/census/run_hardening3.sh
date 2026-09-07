#!/bin/bash
# P3 census v2 hardening — phase 3: matched-DENSITY controls (2026-08-27)
# QSO @ 0.8-1.1 has 92,881 gal; LRG same window/volume has 595,419 (6.4x denser); ELG denser still.
# Subsample LRG and ELG to the QSO count in the SAME window, same scales (15-100):
#   if D2 jumps toward ~2.8 -> sparseness inflates D2 (QSO anomaly explained, amplitude story holds);
#   if D2 stays put        -> QSO census class is genuinely different (real physics-shaped object).
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/hardening-20260827.log
run() { echo "=== census_run.py $* === $(date '+%F %T')" >> "$LOG"; $PY census_run.py "$@" >> "$LOG" 2>&1; }
echo "PHASE-3 START $(date '+%F %T')" >> "$LOG"
run LRG NGC 0.8 1.1 15 100 --nmax 92881 --tag z08-11-nQSO
run ELG_LOPnotqso NGC 0.8 1.1 15 100 --nmax 92881 --tag z08-11-nQSO
echo "PHASE-3 DONE $(date '+%F %T')" >> "$LOG"
