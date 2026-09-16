#!/usr/bin/env python3
"""Resumable H/He/C/N/O correlated atomic controls. See correlated_atoms_notes.md.

Nonrelativistic, infinite-inertia point-component model in calibrated E*,a*.
Spin-pure ROHF reference -> spin-block semicanonical UCCSD(T), all entries.
No measured energy/size/response enters the solver. PySCF 2.14.0.
"""
from __future__ import annotations
import argparse
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
import numpy as np
from scipy.integrate import simpson, cumulative_simpson
import pyscf
from pyscf import cc, fci, gto, lib, scf

REVISION = 'correlated-atom-nr-v1'
Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8}
SPINS = {'H': (1, 0), 'He': (0, 1), 'C': (2, 1), 'N': (3, 2), 'O': (2, 3)}
TERMS = {'H': ('2S','1S'), 'He': ('1S','2S'), 'C': ('3P','2P'),
         'N': ('4S','3P'), 'O': ('3P','4S')}


def dump(path, obj):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basis_for(symbol, level, augmentation):
    if augmentation not in (1, 2, 3):
        raise ValueError('augmentation must be 1, 2, or 3')
    letter = {2:'d', 3:'t', 4:'q'}[level]
    val = gto.basis.load(f'cc-pv{letter}z', symbol)
    aug = gto.basis.load(f'aug-cc-pv{letter}z', symbol)
    assert all(shell in aug for shell in val)
    diffuse = [shell for shell in aug if shell not in val]
    base = gto.basis.load(f'cc-pcv{letter}z', symbol) if Z[symbol]>2 else val
    basis = copy.deepcopy(base+diffuse)
    for extra in range(1, augmentation):
        for shell in diffuse:
            assert len(shell) == 2 and len(shell[1]) == 2
            basis.append([shell[0], [shell[1][0]/(3**extra), 1.]])
    return basis


def occupations(symbol, charge, branch):
    n = Z[symbol]-charge
    if n <= 2:
        return {'A1':((n+1)//2,n//2), 'A2':(0,0), 'B1':(0,0), 'B2':(0,0)}
    # Filled 1s,2s in A1; real p occupations supplied as the term branch.
    p = n-4
    if p == 1: slots = [(0,0),(0,0),(1,0)]
    elif p == 2: slots = [(1,0),(1,0),(0,0)]
    elif p == 3: slots = [(1,0),(1,0),(1,0)]
    elif p == 4: slots = [(1,0),(1,0),(1,1)]
    else: raise ValueError('unsupported sector')
    if branch == 'perp': slots[0],slots[2] = slots[2],slots[0]
    return {'A1':tuple(2+i for i in slots[2]), 'A2':(0,0),
            'B1':slots[0], 'B2':slots[1]}


def radial(mol, density):
    # Exact spherical average of the AO matrix, independent of chosen P axis.
    loc = mol.ao_loc_nr(); labels = np.zeros(mol.nao, dtype=int)
    for b in range(mol.nbas): labels[loc[b]:loc[b+1]]=mol.bas_angular(b)
    avg = np.zeros_like(density)
    for l in np.unique(labels):
        ix=np.flatnonzero(labels==l); g=2*l+1; n=len(ix)//g
        block=density[np.ix_(ix,ix)].reshape(n,g,n,g)
        small=np.einsum('imjm->ij',block)/g
        avg[np.ix_(ix,ix)]=np.einsum('ij,mn->imjn',small,np.eye(g)).reshape(n*g,n*g)
    r=np.geomspace(1e-8,200.,2401); xyz=np.zeros((len(r),3)); xyz[:,2]=r
    ao=mol.eval_gto('GTOval_sph',xyz)
    p=4*np.pi*r*r*np.einsum('pi,ij,pj->p',ao,avg,ao,optimize=True)
    x=np.log(r); count=mol.nelectron
    moments=[float(simpson(p*r**(k+1),x=x)) for k in range(3)]
    cumulative=cumulative_simpson(p*r,x=x,initial=0)
    return {'count':moments[0], 'mean_r':moments[1]/count,
            'mean_r2':moments[2]/count, 'rms_r':float(np.sqrt(moments[2]/count)),
            'r90':float(np.interp(.9*count,cumulative,r)),
            'coarse_mean_difference':abs(float(simpson((p*r*r)[::2],x=x[::2]))-moments[1])/count,
            'r2_integral_check':float(np.einsum('ij,ji',density,mol.intor_symmetric('int1e_r2')))/count,
            'profile':np.column_stack((r[::6],p[::6]/count)).tolist()}


def solve(symbol, level, out, charge=0, branch='parallel', field=0., augmentation=1,
          threads=1, fci_control=False):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    spec={'revision':REVISION,'symbol':symbol,'level':level,'augmentation':augmentation,
          'charge':charge,'branch':branch,'field':field,'fci_control':fci_control}
    key=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    name=f'{symbol}-q{charge}-L{level}-a{augmentation}-{branch}-F{field:+.6f}-{key}'
    record=out/(name+'.json')
    with (out/(name+'.lock')).open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if record.exists():
            saved=json.loads(record.read_text())
            if saved.get('complete') and saved['spec']==spec: return saved
        started=time.time(); lib.num_threads(threads)
        mol=gto.M(atom=f'{symbol} 0 0 0',unit='Bohr',basis=basis_for(symbol,level,augmentation),
                  spin=SPINS[symbol][charge],charge=charge,symmetry='C2v',
                  verbose=4,output=str(out/(name+'.log')),max_memory=1700)
        mf=scf.ROHF(mol)
        mf.irrep_nelec={k:v for k,v in occupations(symbol,charge,branch).items() if k in mol.irrep_name}
        mf.conv_tol=1e-12; mf.conv_tol_grad=1e-8; mf.max_cycle=150
        mf.chkfile=str(out/(name+'.chk'))
        dip=mol.intor_symmetric('int1e_r',comp=3)
        h0=mf.get_hcore(); h=h0+field*dip[2]
        mf.get_hcore=lambda *args:h
        mf.kernel()
        grad=float(np.linalg.norm(mf.get_grad(mf.mo_coeff,mf.mo_occ)))
        if not mf.converged or grad>1e-7:
            raise RuntimeError(f'ROHF not converged: {name}, grad={grad}')
        actual={k:list(v) for k,v in mf.get_irrep_nelec().items()}
        uhf=mf.to_uhf()
        uhf.mo_energy,uhf.mo_coeff=uhf.canonicalize(uhf.mo_coeff,uhf.mo_occ)
        uhf.e_tot=mf.e_tot
        s=mol.intor_symmetric('int1e_ovlp')
        result={'spec':spec,'term_input':TERMS[symbol][charge],'pyscf':pyscf.__version__,
                'source_sha256':digest(Path(__file__)),'basis':mol._basis[symbol],
                'nao':mol.nao,'overlap_min':float(np.linalg.eigvalsh(s)[0]),
                'spin_2Ms':mol.spin,'reference_s2':float(mf.spin_square()[0]),
                'input_occupations':mf.irrep_nelec,'actual_occupations':actual,
                'reference_energy':float(mf.e_tot),'reference_gradient':grad,
                'energies':{'reference':float(mf.e_tot)}}
        arrays={'mo_coeff':np.asarray(uhf.mo_coeff),'mo_occ':uhf.mo_occ,'mo_energy':uhf.mo_energy,
                'overlap':s,'hcore_unperturbed':h0,'dipole':dip}
        if mol.nelectron == 1:
            result['energies'].update(ccsd=float(mf.e_tot),ccsdt=float(mf.e_tot),triples=0.)
            result.update(cc_converged=True,cc_amplitude_residual=0.,t1_norm=0.,t2_max=0.)
            density=uhf.make_rdm1().sum(axis=0)
        else:
            calc=cc.UCCSD(uhf)
            calc.conv_tol=1e-11;calc.conv_tol_normt=1e-9;calc.max_cycle=150
            calc.diis_space=8;calc.max_memory=1700
            eris=calc.ao2mo()
            calc.kernel(eris=eris)
            updated=calc.update_amps(calc.t1,calc.t2,eris)
            residual=float(np.linalg.norm(calc.amplitudes_to_vector(*updated)-calc.amplitudes_to_vector(calc.t1,calc.t2)))
            if not calc.converged or residual>1e-8:
                raise RuntimeError(f'CCSD not converged: {name}, residual={residual}')
            triples=float(calc.ccsd_t(eris=eris)) if mol.nelectron>2 else 0.
            result['energies'].update(ccsd=float(calc.e_tot),ccsdt=float(calc.e_tot+triples),triples=triples)
            result.update(cc_converged=True,cc_amplitude_residual=residual,
                          t1_norm=float(np.sqrt(sum(np.linalg.norm(t)**2 for t in calc.t1))),
                          t1_norm_per_sqrt_entry=float(np.sqrt(sum(np.linalg.norm(t)**2 for t in calc.t1)/mol.nelectron)),
                          t2_max=float(max(np.max(np.abs(t)) for t in calc.t2 if t.size)))
            # Preserve the completed energy before the separate density step.
            dump(out/(name+'.energy.json'),result)
            if field==0:
                calc.solve_lambda(eris=eris)
                if not calc.converged_lambda: raise RuntimeError('CCSD lambda did not converge')
                d=calc.make_rdm1(ao_repr=True)
                density=d[0]+d[1]
                for family,parts in [('t1',calc.t1),('t2',calc.t2),('l1',calc.l1),('l2',calc.l2)]:
                    for i,part in enumerate(parts):arrays[f'{family}_{i}']=part
                result['lambda_converged']=bool(calc.converged_lambda)
            if fci_control and mol.nelectron==2 and field==0:
                from pyscf import ao2mo
                c=mf.mo_coeff
                fcisolver=fci.direct_spin1.FCI(mol);fcisolver.conv_tol=1e-12
                efci,_=fcisolver.kernel(c.T@h@c,ao2mo.kernel(mol,c),mol.nao,mol.nelec)
                result['fci_energy']=float(efci)
                result['ccsd_minus_fci']=float(calc.e_tot-efci)
        if field==0:
            result['density_count']=float(np.einsum('ij,ji',density,s))
            result['radial_ccsd']=radial(mol,density)
            arrays['density_ccsd']=density
            state=out/(name+'.npz');np.savez_compressed(state,**arrays)
            result['state_file']=state.name;result['state_sha256']=digest(state)
        result['elapsed_seconds']=time.time()-started
        result['complete']=True
        dump(record,result)
        print(json.dumps({'saved':record.name,'energy':result['energies']['ccsdt'],
                          'seconds':round(result['elapsed_seconds'],2)}),flush=True)
        mol.stdout.close()
        return result


def campaign(args):
    # Fixed order: every basis has neutral/removal energies before its responses.
    branches=['parallel','perp'] if args.symbol in ('C','O') else ['parallel']
    for level in args.levels:
        for augmentation in args.augmentations:
            solve(args.symbol,level,args.out,threads=args.threads,augmentation=augmentation,
                  fci_control=args.symbol=='He')
            if args.symbol in ('C','N','O','He'):
                solve(args.symbol,level,args.out,charge=1,threads=args.threads,augmentation=augmentation)
            if not args.energies_only:
                for branch in branches:
                    if branch!='parallel':
                        solve(args.symbol,level,args.out,branch=branch,threads=args.threads,augmentation=augmentation)
                    for field in (0.001,-0.001,0.0005,-0.0005):
                        solve(args.symbol,level,args.out,branch=branch,field=field,
                              threads=args.threads,augmentation=augmentation)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--symbol',choices=Z,required=True)
    p.add_argument('--levels',type=int,nargs='+',default=[2,3,4])
    p.add_argument('--augmentations',type=int,nargs='+',default=[1])
    p.add_argument('--out',required=True)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--energies-only',action='store_true')
    campaign(p.parse_args())
