#!/bin/bash
# P3 census v2 hardening — phase 2: multi-random on the big tracers (2026-08-27)
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/hardening-20260827.log
run() { echo "=== census_run.py $* === $(date '+%F %T')" >> "$LOG"; $PY census_run.py "$@" >> "$LOG" 2>&1; }
echo "PHASE-2 START $(date '+%F %T')" >> "$LOG"
run LRG NGC 0.4 1.1 10 100 --ran 1 --tag ran1
run ELG_LOPnotqso NGC 0.8 1.6 10 100 --ran 1 --tag ran1
echo "PHASE-2 DONE $(date '+%F %T')" >> "$LOG"
