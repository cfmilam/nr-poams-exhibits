#!/usr/bin/env python3
"""Run independent correlated-atom jobs with restartable per-field records."""
import argparse
import concurrent.futures
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--symbols',nargs='+',default=['H','He','C','N','O'])
    p.add_argument('--levels',nargs='+',default=['2','3','4'])
    p.add_argument('--augmentations',nargs='+',default=['1'])
    p.add_argument('--workers',type=int,default=3)
    p.add_argument('--energies-only',action='store_true')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    lock=(a.out/'campaign.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state={'pid':os.getpid(),'started':time.time(),'status':'running','jobs':{}}
    guard=threading.Lock()
    def save():
        tmp=a.out/'campaign-status.tmp';tmp.write_text(json.dumps(state,indent=2)+'\n')
        tmp.replace(a.out/'campaign-status.json')
    save()
    def work(symbol):
        cmd=[sys.executable,str(Path(__file__).with_name('correlated_atoms.py')),
             '--symbol',symbol,'--out',str(a.out/'records'),'--levels',*a.levels,
             '--augmentations',*a.augmentations]
        if a.energies_only:cmd.append('--energies-only')
        with (a.out/(symbol+'-worker.log')).open('a') as log:
            child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
            with guard:
                state['jobs'][symbol]={'pid':child.pid,'started':time.time(),'status':'running'};save()
            code=child.wait()
            with guard:
                state['jobs'][symbol].update(exit_code=code,finished=time.time(),
                                            status='completed' if code==0 else 'failed');save()
            return code
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        codes=list(pool.map(work,a.symbols))
    state['status']='completed' if all(c==0 for c in codes) else 'needs_attention'
    state['finished']=time.time();save()


if __name__=='__main__':main()
