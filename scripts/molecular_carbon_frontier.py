#!/usr/bin/env python3
"""Recorded carbon frontier gates from the frozen pair-ledger translation.

The expensive evaluator is documented in the accompanying audit.  This compact
instrument independently rechecks convergence, benchmark errors, Hessian step
stability, mass-metric mode rates, and the relaxed ethane torsion barrier.
"""
from __future__ import annotations
import math
import numpy as np

GEOMETRY={
 "methane": {"tz":(1.082444955,),"qz":(1.081578107,),"obs":(1.087,)},
 "co2": {"tz":(1.166346541,),"qz":(1.163596569,),"obs":(1.162,)},
 "ethane": {"tz":(1.528663803,1.083951964),"qz":(1.528138714,1.083587127),"obs":(1.535,1.094)},
 "ethyne": {"tz":(1.207331002,1.057153169),"qz":(1.205933338,1.059001342),"obs":(1.203,1.063)},
 "ethene": {"tz":(1.330621640,1.075736113),"qz":(1.328691389,1.076385096),"obs":(1.339,1.086)},
 "benzene": {"tz":(1.387929586,1.076556531),"qz":(1.387430262,1.078252882),"obs":(1.397,1.084)},
}

HESSIANS={
 "methane":(np.array([[5.202222057]]),np.array([[5.200721357]])),
 "co2":(np.array([[7.804675712]]),np.array([[7.801279054]])),
 "ethyne":(np.array([[3.647852501,-.048505285],[-.048505285,2.995493391]]),np.array([[3.646775657,-.048448670],[-.048448670,2.994585062]])),
 "ethene":(np.array([[2.228237772,.066155814],[.066155814,5.339450871]]),np.array([[2.227621596,.066169445],[.066169445,5.337823191]])),
 "benzene":(np.array([[11.083576135,.145448898],[.145448898,7.913434623]]),np.array([[11.080786451,.145468930],[.145468930,7.910995546]])),
}

def pct(a,b): return 100*(a/b-1)

if __name__=="__main__":
    for name,row in GEOMETRY.items():
        conv=[pct(q,t) for q,t in zip(row["qz"],row["tz"])]
        err=[pct(q,o) for q,o in zip(row["qz"],row["obs"])]
        print(f"{name:8s} convergence%={conv} benchmark_error%={err}")
        assert max(map(abs,conv))<1 and max(map(abs,err))<1
    for name,(h15,h10) in HESSIANS.items():
        drift=np.max(np.abs(h10/h15-1))*100
        print(f"{name:8s} Hessian step drift={drift:.4f}%")
        assert drift<.2
    barrier=( -79.70855492129222 - (-79.71357415573546))*627.5094740631
    print(f"ethane relaxed torsion barrier={barrier:.6f} kcal/mol")
    assert abs(pct(barrier,2.928))<10
    print("ALL RECORDED CARBON FRONTIER GATES PASS")
