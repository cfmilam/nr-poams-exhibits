#!/bin/bash
# P3 census hardening round 5 (2026-08-29) — QSO anomaly next-tier: imaging systematics + fiber assignment
# Cuts (executable tonight, no imaging maps needed):
#   (1) QSO cross-hemisphere: SGC has INDEPENDENT imaging vs NGC.
#       If QSO_SGC D2 ~ QSO_NGC D2 (2.795) => imaging-systematics patterns DISFAVORED.
#       Large split => points AT imaging.
#   (2) QSO within-class z-ladder: 1.1-1.45 and 1.45-1.8 (NGC). Imaging systematics scale
#       with depth/target density; FLAT ladder further disfavors them.
#   (3) LRG SGC matched-z reference: gives SGC its own LRG<->QSO comparison pair.
cd "$(dirname "$0")"
PY=venv/bin/python
LOG=results/hardening5-20260829.log
run() { echo "=== census_run.py $* === $(date '+%F %T')" >> "$LOG"; $PY census_run.py "$@" >> "$LOG" 2>&1; }

echo "HARDENING-5 DRIVER START $(date '+%F %T')" >> "$LOG"

# (1) QSO cross-hemisphere — SGC, matched z-window/scales to NGC baseline
run QSO SGC 0.8 1.1 15 100 --tag z08-11

# (2) QSO within-class z-ladder (NGC), fit 15-100
run QSO NGC 1.1 1.45 15 100 --tag z11-145
run QSO NGC 1.45 1.8 15 100 --tag z145-18

# (3) LRG SGC matched-z reference (data already in data/)
run LRG SGC 0.8 1.1 15 100 --tag z08-11

echo "HARDENING-5 DRIVER DONE $(date '+%F %T')" >> "$LOG"
echo "ALL DONE" >> "$LOG"
