#!/bin/bash
# P3 census v2 hardening — sequential driver (2026-08-27)
# Items from census-run-results.md v2 list executable tonight without new downloads:
#   (1) matched-z cross-tracer trio (LRG/ELG/QSO at z 0.8-1.1, common fit range 15-100)
#       -> separates depth-vs-tracer-bias readings of the D2 ladder
#   (2) within-tracer z ladders (LRG 3 slices, ELG 2, BGS 2)
#       -> if D2 rises WITH z inside one tracer, depth is doing work; if flat inside but
#          jumping across tracers, bias is doing work
#   (3) estimator systematics on BGS NGC: second random realization (--ran 1),
#       embedding-criterion variants (--fill 0.90 / 0.99)
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/hardening-20260827.log
run() { echo "=== census_run.py $* === $(date '+%F %T')" >> "$LOG"; $PY census_run.py "$@" >> "$LOG" 2>&1; }

echo "HARDENING DRIVER START $(date '+%F %T')" >> "$LOG"

# --- (3) cheap estimator systematics first (BGS is the fast tracer) ---
run BGS_BRIGHT-21.5 NGC 0.05 0.40 5 70 --ran 1 --tag ran1
run BGS_BRIGHT-21.5 NGC 0.05 0.40 5 70 --fill 0.90 --tag fill090
run BGS_BRIGHT-21.5 NGC 0.05 0.40 5 70 --fill 0.99 --tag fill099

# --- (1) matched-z cross-tracer trio (the bias discriminator) ---
run LRG NGC 0.8 1.1 15 100 --tag z08-11
run ELG_LOPnotqso NGC 0.8 1.1 15 100 --tag z08-11
run QSO NGC 0.8 1.1 15 100 --tag z08-11

# --- (2) within-tracer ladders ---
run BGS_BRIGHT-21.5 NGC 0.05 0.22 5 70 --tag z005-022
run BGS_BRIGHT-21.5 NGC 0.22 0.40 5 70 --tag z022-040
run LRG NGC 0.4 0.6 10 100 --tag z04-06
run LRG NGC 0.6 0.8 10 100 --tag z06-08
run ELG_LOPnotqso NGC 1.1 1.6 10 100 --tag z11-16

echo "ALL DONE $(date '+%F %T')" >> "$LOG"
