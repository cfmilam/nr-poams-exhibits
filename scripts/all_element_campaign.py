#!/usr/bin/env python3
"""Checkpointed two-level campaign; independent processes, not model agents."""
import argparse
import concurrent.futures
import json
import os
import time
import traceback
from pathlib import Path

os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')

def job(z,level,out):
    # A manager-specific marker drains queued work at job boundaries. Active
    # atomic/removal calculations finish normally before a pool is resized.
    # Remove the marker after that manager exits; never leave it for PID reuse.
    if Path(out,f'campaign-pause-{os.getppid()}.json').is_file():
        return {'z':z,'level':level,'completed':False,'paused':True}
    from all_element_atoms import calculate, removal_calculation
    try:
        calculate(z,level,out)
        if level==3:
            matches=list(Path(out).glob(f'{z:03d}-*-v3z.json'))
            if matches:removal_calculation(matches[0])
        return {'z':z,'level':level,'completed':True}
    except Exception:
        result={'z':z,'level':level,'completed':False,'error':traceback.format_exc()}
        Path(out,f'{z:03d}-v{level}z-error.json').write_text(json.dumps(result,indent=2))
        return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--start',type=int,default=1)
    parser.add_argument('--end',type=int,default=118)
    parser.add_argument('--levels',type=int,nargs='+',default=[2,3])
    args=parser.parse_args()
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    order=[z for z in [1,2,6,7,8]+list(range(3,119)) if args.start<=z<=args.end]
    order=list(dict.fromkeys(order))
    begun=time.time();ledger=[]
    for level in args.levels:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,max_tasks_per_child=1) as pool:
            futures=[pool.submit(job,z,level,str(out)) for z in order]
            for future in concurrent.futures.as_completed(futures):
                entry=future.result();ledger.append(entry)
                Path(out,f'campaign-progress-{os.getpid()}.json').write_text(json.dumps({'started':begun,'updated':time.time(),'jobs':ledger},indent=2))
                print('CAMPAIGN',json.dumps({**entry,'elapsed_minutes':round((time.time()-begun)/60,2)}),flush=True)
