#!/usr/bin/env python3
"""Independent saved-result identities and explicitly graded numerical checks."""
import argparse
import ast
import json
import math
import os
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import numpy as np
from pyscf import gto

def verify(source,require_all=False):
    source=Path(source);checks=[];issues=[];states=[];reduced_ranks=[]
    def check(condition,label):
        if not condition:raise AssertionError(label)
        checks.append(label)
    files=sorted(source.glob('*-v?z.json'))
    for file in files:
        r=json.loads(file.read_text());states.append(r);tag=f"{r['symbol']}/v{r['level']}z"
        if not r['validated']:
            issues.append({'state':tag,'issue':'state not validated'});continue
        check(abs(r['count']-r['z'])<2e-7,tag+' entry conservation')
        check(len(r['occupied_orbitals'])==r['z'],tag+' occupied entries')
        ns=[sum(o['spin']==s for o in r['occupied_orbitals']) for s in [0,1]]
        check(ns[0]-ns[1]==r['input_spin_2Ms'],tag+' specified spin projection')
        check(r['orthonormality_error']<2e-7,tag+' orbital metric')
        check(r['gradient_norm']<3e-5,tag+' stationary RI-JK mean-field state')
        check(math.isfinite(r['direct_gradient_norm']) and math.isfinite(r['integral_energy_difference']),tag+' finite direct-integral diagnostics')
        if r['direct_gradient_norm']>=3e-5:
            issues.append({'state':tag,'issue':'direct-integral gradient exceeds the RI-JK stationarity threshold',
                           'gradient':r['direct_gradient_norm'],'energy_difference':r['integral_energy_difference']})
        statefile=file.with_name(file.stem+'-state.npz')
        check(statefile.is_file(),tag+' reusable state present')
        with np.load(statefile,allow_pickle=False) as saved:
            coeff=saved['mo_coeff'];energy=saved['mo_energy'];occupation=saved['mo_occ']
            n=r['nao']
            # PySCF 2.14 removes overlap eigenvectors <= 1e-6 by default.
            # The AO dimension and retained orbital dimension need not agree.
            nmo=energy.shape[1] if energy.ndim==2 else 0
            check(0<nmo<=n and coeff.shape==(2,n,nmo) and
                  energy.shape==occupation.shape==(2,nmo),tag+' saved-state dimensions')
            check(all(np.all(np.isfinite(x)) for x in [coeff,energy,occupation]),tag+' finite saved state')
            check(np.all((occupation==0)|(occupation==1)) and int(occupation.sum())==r['z'],tag+' saved-state entry count')
            check(int(occupation[0].sum()-occupation[1].sum())==r['input_spin_2Ms'],tag+' saved-state spin projection')
            check(all(occupation[o['spin'],o['index']]==1 and
                      abs(energy[o['spin'],o['index']]-o['energy'])<1e-12
                      for o in r['occupied_orbitals']),tag+' exported orbitals match reference')
            metadata=json.loads(str(saved['mol_json']));basis=json.loads(str(saved['basis_json']))
            # PySCF dumps instance attributes; singlets may inherit Mole.spin=0.
            check(metadata['charge']==0 and metadata.get('spin',0)==r['input_spin_2Ms'] and
                  ast.literal_eval(metadata['basis'])==basis and
                  ast.literal_eval(metadata['atom']).split()[0]==r['symbol'],
                  tag+' saved neutral input metadata')
            # Reconstruct the AO metric from the exported molecule, rather than
            # trusting the scalar normalization diagnostics in the source JSON.
            molecule=gto.loads(str(saved['mol_json']))
            overlap=molecule.intor_symmetric('int1e_ovlp')
            overlap_values,overlap_vectors=np.linalg.eigh(overlap)
            keep=overlap_values>1e-6
            check(nmo==int(keep.sum()),tag+' retained dimension matches overlap cutoff')
            restored_metric=max(float(np.max(np.abs(c.T@overlap@c-np.eye(nmo)))) for c in coeff)
            if nmo<n:
                discarded_projection=max(float(np.max(np.abs(overlap_vectors[:,~keep].T@c))) for c in coeff)
                check(discarded_projection<2e-7,tag+' no discarded overlap-mode coefficients')
                record={'state':tag,'nao':n,'nmo_per_spin':nmo,'overlap_cutoff':1e-6,
                        'discarded_overlap_eigenvalues':overlap_values[~keep].tolist(),
                        'first_retained_overlap_eigenvalue':float(overlap_values[keep][0]),
                        'maximum_discarded_mode_coefficient':discarded_projection}
                reduced_ranks.append(record)
                issues.append(dict(record,issue='linearly dependent AO modes removed by eigensolver'))
            restored_count=sum(float(np.einsum('ij,ji',(c*o)@c.T,overlap)) for c,o in zip(coeff,occupation))
            check(restored_metric<2e-7,tag+' independently restored AO metric')
            check(abs(restored_count-r['z'])<2e-7,tag+' independently restored density count')
        check(abs(r['one_body_energy']+r['direct_energy']+r['exchange_energy']-r['energy'])<2e-7,tag+' energy account')
        radial=r['radial']
        check(abs(radial['normalization']-r['z'])<2e-5,tag+' independent radial entry integral')
        check(abs(radial['matrix_normalization']-r['z'])<2e-7,tag+' X2C metric reconstruction')
        check(abs(radial['normalization']-radial['coarse_normalization'])<2e-5,tag+' radial grid convergence')
        check(0<radial['mean_r']<=radial['rms_r'] and radial['r90']>0,tag+' radius moment inequality')
        resp=r['response']
        if 'frozen_tensor' in resp:
            check(abs(np.trace(resp['frozen_tensor'])/3-resp['frozen_spectral_static'])<2e-7,
                  tag+' frozen spectrum versus static tensor')
        if resp.get('certified'):
            t=np.array(resp['relaxed_tensor'])
            check(np.max(np.abs(t-t.T))<2e-5 and np.min(np.linalg.eigvalsh(t))>0,tag+' positive reciprocal response')
            check(resp['scaled_linear_residual']<2e-8,tag+' independent response equation residual')
            check(resp['tensor_asymmetry']<2e-5,tag+' reciprocity before tensor symmetrization')
        else:
            issues.append({'state':tag,'issue':resp['status'],'response':resp.get('relaxed_isotropic')})
        expected=r['input_configuration_l_counts']
        actual=[sum(o['l']==l for o in r['occupied_orbitals']) for l in range(4)]
        if expected!=actual:issues.append({'state':tag,'issue':'angular configuration changed during relaxation','input':expected,'output':actual})
    pairs={(r['z'],r['level']) for r in states if r['validated']}
    if require_all:
        check(all((z,l) in pairs for z in range(1,119) for l in [2,3]),'118 validated states at both basis sizes')
    ions=[];deferred_ions=[]
    neutral_by_pair={(r['z'],r['level']):r for r in states}
    for file in sorted(source.glob('*-removal.json')):
        r=json.loads(file.read_text());pair=(r['z'],r['level'])
        # Workers can finish another atom after the neutral file list was read.
        # Keep this audit a fixed neutral snapshot, reporting newer ions separately.
        if pair not in neutral_by_pair:
            deferred_ions.append(file.name);continue
        ions.append(r);neutral=neutral_by_pair[pair];tag=f"Z={r['z']}/v{r['level']}z"
        attempts=r['attempts']
        if r['z']==1:
            check(not attempts and r['certified'] and r['ion_energy']==0 and
                  r['ion_spin_2Ms']==0,tag+' exact empty-ion control')
        else:
            ns=[sum(o['spin']==spin for o in neutral['occupied_orbitals']) for spin in [0,1]]
            expected_seeds={spin for spin in [0,1] if ns[spin]>0}
            check(len(attempts)==len(expected_seeds) and
                  {a['removed_spin'] for a in attempts}==expected_seeds,
                  tag+' every occupied spin-removal seed recorded once')
            for a in attempts:
                spin=a['removed_spin'];counts=ns.copy();counts[spin]-=1
                check(a['ion_spin_2Ms']==abs(counts[0]-counts[1]) and
                      a['ion_spin_2Ms']<=r['z']-1 and
                      (r['z']-1-a['ion_spin_2Ms'])%2==0,
                      tag+f' removal seed {spin} entry/spin bookkeeping')
                check(all(math.isfinite(a[k]) for k in ['energy','gradient_norm','spin_squared']) and
                      a['gradient_norm']>=0,tag+f' removal seed {spin} finite diagnostics')
                check(bool(a['history']) and abs(a['history'][-1]['energy']-a['energy'])<1e-10,
                      tag+f' removal seed {spin} solver history energy')
            eligible=[a for a in attempts if a['converged'] and a['gradient_norm']<3e-5]
            check(bool(eligible)==r['certified'],tag+' removal certification follows seed stationarity')
            if len(eligible)<len(attempts):
                issues.append({'state':tag,'issue':'not every ionic seed reached stationarity',
                               'eligible':len(eligible),'attempted':len(attempts)})
            if eligible:
                minimum=min(a['energy'] for a in eligible)
                check(abs(r['ion_energy']-minimum)<1e-10 and
                      any(abs(a['energy']-minimum)<1e-10 and a['ion_spin_2Ms']==r['ion_spin_2Ms'] for a in eligible),
                      tag+' selected ion is lowest stationary tested seed')
        if not r['certified']:
            issues.append({'state':tag,'issue':'removal interval not certified'});continue
        check(abs(r['ion_energy']-neutral['energy']-r['interval'])<1e-10,tag+' separate-state removal account')
        check(r['interval']>0,tag+' bound neutral against tested removal')
    if require_all:
        check({r['z'] for r in ions if r['level']==3}==set(range(1,119)),
              '118 larger-basis removal calculations attempted and recorded')
    fine=[r for r in states if r['level']==3 and r['validated']]
    basis=[]
    for hi in fine:
        lo=next((r for r in states if r['z']==hi['z'] and r['level']==2 and r['validated']),None)
        if not lo:continue
        diff=None
        if hi['response'].get('certified') and lo['response'].get('certified'):
            diff=100*(hi['response']['relaxed_isotropic']/lo['response']['relaxed_isotropic']-1)
        basis.append({'z':hi['z'],'symbol':hi['symbol'],'energy_change':hi['energy']-lo['energy'],
                      'mean_radius_percent':100*(hi['radial']['mean_r']/lo['radial']['mean_r']-1),
                      'response_percent':diff})
        if diff is not None and abs(diff)>5:issues.append({'state':hi['symbol'],'issue':'response basis sensitivity exceeds 5%','percent':diff})
    result={'checks_passed':len(checks),'checks':checks,'state_results':len(states),'fine_states':len(fine),
            'fine_responses':sum(r['response'].get('certified',False) for r in fine),
            'fine_removal_records':sum(r['level']==3 for r in ions),
            'fine_removal_intervals':sum(r['level']==3 and r['certified'] for r in ions),
            'deferred_removal_records':deferred_ions,
            'basis_comparisons':basis,'reduced_rank_states':reduced_ranks,'issues':issues,
            'max_direct_gradient_norm':max(r['direct_gradient_norm'] for r in states),
            'max_direct_integral_energy_error':max(abs(r['integral_energy_difference']) for r in states),
            'warning':'Passed numerical checks do not establish correlated accuracy, physical ground-term selection or internal core structure.'}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--out',required=True);p.add_argument('--require-all',action='store_true');a=p.parse_args()
    r=verify(a.source,a.require_all);Path(a.out).write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k not in ['checks','basis_comparisons']}))
