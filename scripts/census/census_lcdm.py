#!/usr/bin/env python3
"""lcdm-only wrapper around census_run.main (memory/time-lean for the gamma grid)."""
import sys, census_run
census_run.DMAPS = {"lcdm": census_run.d_lcdm}  # restrict to one map
if __name__ == "__main__":
    argv = sys.argv[1:]
    opts, skip = {}, set()
    for i, a in enumerate(argv):
        if a == "--tag": opts["tag"] = argv[i+1]; skip.update({i, i+1})
    args = [a for i, a in enumerate(argv) if i not in skip]
    tracer, cap = args[0], args[1]
    zmin, zmax, rmin, rmax = map(float, args[2:6])
    census_run.main(tracer, cap, zmin, zmax, rmin, rmax, **opts)
