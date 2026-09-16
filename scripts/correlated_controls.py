#!/usr/bin/env python3
"""Independent all-entry FCI vs UCCSD(T) in small cc-pVDZ C/N/O spaces.

This measures correlation truncation in that finite space only. It is not a
basis-limit error estimate. No frozen entries or fitted target energies.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from pyscf import ao2mo,cc,fci,gto,lib,scf
from correlated_atoms import occupations,SPINS,dump,digest


def run(symbol,charge,out):
    path=out/f'{symbol}-q{charge}-cc-pvdz-fci.json'
    if path.exists():return
    start=time.time();lib.num_threads(1)
    mol=gto.M(atom=f'{symbol} 0 0 0',basis='cc-pvdz',spin=SPINS[symbol][charge],
              charge=charge,symmetry='C2v',verbose=4,output=str(path.with_suffix('.log')),max_memory=1600)
    mf=scf.ROHF(mol);mf.irrep_nelec=occupations(symbol,charge,'parallel')
    mf.conv_tol=1e-12;mf.conv_tol_grad=1e-9;mf.kernel()
    assert mf.converged
    uhf=mf.to_uhf();uhf.mo_energy,uhf.mo_coeff=uhf.canonicalize(uhf.mo_coeff,uhf.mo_occ)
    c=cc.UCCSD(uhf);c.conv_tol=1e-11;c.conv_tol_normt=1e-9;c.max_cycle=150
    eris=c.ao2mo();c.kernel(eris=eris);assert c.converged
    triples=float(c.ccsd_t(eris=eris))
    ci=fci.direct_spin1.FCI(mol);ci.conv_tol=1e-11;ci.max_cycle=180;ci.max_space=24
    coeff=mf.mo_coeff
    e,v=ci.kernel(coeff.T@mf.get_hcore()@coeff,ao2mo.kernel(mol,coeff),mol.nao,mol.nelec)
    assert ci.converged
    spin2=float(ci.spin_square(v,mol.nao,mol.nelec)[0])
    result={'symbol':symbol,'charge':charge,'model':'nonrelativistic point compact; infinite inertia',
            'basis':mol._basis[symbol],'basis_name':'cc-pVDZ','entries':mol.nelectron,
            'orbitals':mol.nao,'fci_converged':bool(ci.converged),'cc_converged':bool(c.converged),
            'reference_energy':float(mf.e_tot),'ccsd_energy':float(c.e_tot),
            'ccsdt_energy':float(c.e_tot+triples),'fci_energy':float(e),'triples':triples,
            'ccsdt_minus_fci':float(c.e_tot+triples-e),'fci_s2':spin2,
            'source_sha256':digest(Path(__file__)),'seconds':time.time()-start}
    dump(path,result);print(json.dumps(result,default=str),flush=True)
    mol.stdout.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    for symbol in ('C','N','O'):
        for charge in (0,1):run(symbol,charge,a.out)
