#!/bin/bash
# P3 census hardening round 4 (2026-08-28) — QSO anomaly kill-tests
# From census-run-results.md "Mundane suspects remaining": catastrophic-z, weight handling.
#   (1) QSO weight-handling check: --noweight rerun at matched-z window
#   (2) Catastrophic-z injection response curve on LRG (the clean class):
#       scramble f in {2%, 5%, 10%, 20%} at matched-z 0.8-1.1, fit 15-100.
#       Question: what f is needed to fake QSO's +0.27 (15-100) / +0.20 (30-100)?
#       DESI QSO catastrophic-z rate ~1-3% => if needed f >> 3%, suspect is DEAD.
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/hardening4-20260828.log
run() { echo "=== census_run.py $* === $(date '+%F %T')" >> "$LOG"; $PY census_run.py "$@" >> "$LOG" 2>&1; }

echo "HARDENING-4 DRIVER START $(date '+%F %T')" >> "$LOG"

# (1) weight handling — QSO is the fast tracer (93k), do first
run QSO NGC 0.8 1.1 15 100 --noweight --tag z08-11-nowt

# (2) catastrophic-z response curve on LRG at the matched-z window
run LRG NGC 0.8 1.1 15 100 --scramble 0.02 --tag z08-11-scr02
run LRG NGC 0.8 1.1 15 100 --scramble 0.05 --tag z08-11-scr05
run LRG NGC 0.8 1.1 15 100 --scramble 0.10 --tag z08-11-scr10
run LRG NGC 0.8 1.1 15 100 --scramble 0.20 --tag z08-11-scr20

# (2b) sanity anchor: scramble applied to QSO itself at its own catastrophic rate —
#      if QSO's D2 barely moves at 2%, in-window scrambling can't be its own story either
run QSO NGC 0.8 1.1 15 100 --scramble 0.02 --tag z08-11-scr02

echo "HARDENING-4 DRIVER DONE $(date '+%F %T')" >> "$LOG"
echo "ALL DONE" >> "$LOG"
