#!/usr/bin/env python3
"""Audit saved correlated atomic calculations and build the exhibit handoff."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
from correlated_atoms import SPINS,Z,dump

EXPERIMENT={'C':11.2602880,'N':14.53413,'O':13.618055}
EXPERIMENT_UNC={'C':.0000011,'N':.00004,'O':.000007}
THEORY={'C':11.67,'N':7.26,'O':5.24}
THEORY_COMPONENTS={'C':[10.27,12.38],'N':[7.26,7.26],'O':[5.86,4.94]}
NIST='https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&format=1&level_out=on&seq_out=on&shells_out=on&unc_out=on&units=1'
PAPER='https://doi.org/10.1088/0953-4075/31/10/011'


def build(roots,controls,out,strict=False):
    records={}; paths=[]; statepaths=[]; errors=[]
    def check(condition,message):
        if not condition:errors.append(message)
    for root in roots:
        for path in sorted(root.glob('*.json')):
            if '.energy.' in path.name:continue
            r=json.loads(path.read_text())
            if not r.get('complete'):
                check(not strict,'incomplete record '+path.name)
                continue
            s=r['spec'];key=(s['symbol'],s['level'],s['augmentation'],s['charge'],s['branch'],s['field'])
            check(key not in records,f'duplicate specification {key}');records[key]=r;paths.append(path)
            tag=path.name
            check(r['reference_gradient']<1e-7,'reference gradient '+tag)
            check(r['cc_converged'] and r['cc_amplitude_residual']<1e-8,'CC residual '+tag)
            spin=SPINS[s['symbol']][s['charge']]/2
            check(abs(r['reference_s2']-spin*(spin+1))<1e-8,'reference spin '+tag)
            check(r['input_occupations']==r['actual_occupations'],'irrep occupation '+tag)
            check(r['overlap_min']>1e-8,'overlap conditioning '+tag)
            if 'ccsd_minus_fci' in r:check(abs(r['ccsd_minus_fci'])<1e-8,'He FCI '+tag)
            if s['field']==0:
                n=Z[s['symbol']]-s['charge'];rad=r['radial_ccsd']
                check(abs(r['density_count']-n)<1e-7,'density normalization '+tag)
                check(abs(rad['count']-n)<1e-6,'radial normalization '+tag)
                check(abs(rad['mean_r2']-rad['r2_integral_check'])<1e-7,'radial analytic moment '+tag)
                check(rad['coarse_mean_difference']<1e-7,'radial quadrature '+tag)
                sp=path.parent/r['state_file'];statepaths.append(sp)
                check(hashlib.sha256(sp.read_bytes()).hexdigest()==r['state_sha256'],'state hash '+tag)
                with np.load(sp,allow_pickle=False) as state:
                    for c in state['mo_coeff']:
                        check(np.max(np.abs(c.T@state['overlap']@c-np.eye(c.shape[1])))<1e-8,'orbital metric '+tag)
                    check(np.max(np.abs(state['density_ccsd']-state['density_ccsd'].T))<1e-8,'density symmetry '+tag)
    source=Path(__file__).parent
    source_versions={hashlib.sha256((source/'correlated_atoms.py').read_bytes()).hexdigest():source/'correlated_atoms.py'}
    for path in sorted((source/'correlated_atoms_sources').glob('*.py')):
        check(hashlib.sha256(path.read_bytes()).hexdigest()==path.stem,'archived source hash '+path.name)
        source_versions[path.stem]=path
    for r in records.values():
        check(r['source_sha256'] in source_versions,'missing solver source '+r['source_sha256'])
    if strict:
        groups={key[:3] for key in records}
        required={(a,l,1) for a in Z for l in (2,3,4)}
        required.update((a,l,2) for a in Z for l in (3,4))
        for group in sorted(required-groups):check(False,'missing planned basis group '+str(group))
        for a,l,u in sorted(groups):
            branches=('parallel','perp') if a in ('C','O') else ('parallel',)
            expected={(a,l,u,0,b,f) for b in branches for f in (0.,.001,-.001,.0005,-.0005)}
            if a!='H':expected.add((a,l,u,1,'parallel',0.))
            for key in sorted(expected-records.keys()):check(False,'missing planned state/field '+str(key))
    original_controls=json.loads((source/'whole_atom_handoff.json').read_text())
    units=original_controls['units']
    rate=units['rate_Hz']/1e15;length=units['length_m']*1e12
    def get(a,l,u,q=0,b='parallel',f=0.):return records.get((a,l,u,q,b,f))
    def response(a,l,u,b,method):
        rr=[get(a,l,u,b=b,f=f) for f in (0,.001,-.001,.0005,-.0005)]
        if not all(rr):return None
        e=[r['energies'][method] for r in rr]
        large=-(e[1]+e[2]-2*e[0])/1e-6
        small=-(e[3]+e[4]-2*e[0])/2.5e-7
        extrap=(4*small-large)/3
        even=max(abs(e[1]-e[2]),abs(e[3]-e[4]))
        step=abs(small-large)/abs(extrap)
        check(extrap>0,f'positive response {a,l,u,b,method}')
        check(even<1e-8,f'field parity {a,l,u,b,method}')
        check(step<.002,f'field step {a,l,u,b,method}')
        return {'value':extrap,'field_step_relative':step,'signed_energy_asymmetry':even,
                'coarse_curvature':large,'fine_curvature':small}
    atoms={}
    for a in Z:
        levels=[]
        for l in (2,3,4):
            for u in (1,2,3):
                z=get(a,l,u)
                if not z:continue
                ion=get(a,l,u,q=1)
                row={'level':l,'augmentation':u,'basis':{1:'aug',2:'double-aug',3:'triple-aug'}[u]+f'-cc-p{"CV" if Z[a]>2 else "V"}{ {2:"D",3:"T",4:"Q"}[l]}Z',
                     'nao':z['nao'],'energies':z['energies'],'term':z['term_input'],
                     'mean_pm':z['radial_ccsd']['mean_r']*length,
                     'rms_pm':z['radial_ccsd']['rms_r']*length,'r90_pm':z['radial_ccsd']['r90']*length,
                     'state_file':z['state_file'],'state_sha256':z['state_sha256'],
                     'radial_profile':[[r*length,p/length] for r,p in z['radial_ccsd']['profile']],
                     'cc_residual':z['cc_amplitude_residual'],'t1_norm_per_sqrt_entry':z.get('t1_norm_per_sqrt_entry',0),
                     'removal_eV':None,'removal_rate_PHz':None,'response':{}}
                if ion:
                    delta=ion['energies']['ccsdt']-z['energies']['ccsdt']
                    row.update(removal_eV=delta*units['energy_eV'],removal_rate_PHz=delta*rate,
                               reference_removal_eV=(ion['energies']['reference']-z['energies']['reference'])*units['energy_eV'],
                               ccsd_removal_eV=(ion['energies']['ccsd']-z['energies']['ccsd'])*units['energy_eV'])
                    check(delta>0,f'positive removal {a,l,u}')
                for method in ('reference','ccsd','ccsdt'):
                    axial=response(a,l,u,'parallel',method)
                    other=response(a,l,u,'perp',method) if a in ('C','O') else axial
                    if axial and other:
                        row['response'][method]={'parallel':axial,'perpendicular':other,
                                                'average':(axial['value']+2*other['value'])/3}
                perpendicular=get(a,l,u,b='perp')
                if perpendicular:
                    gap=abs(z['energies']['ccsdt']-perpendicular['energies']['ccsdt'])
                    check(gap<1e-8,f'rotated zero-field degeneracy {a,l,u}')
                    row['zero_field_rotation_difference']=gap
                levels.append(row)
        if not levels:continue
        complete=[v for v in levels if 'ccsdt' in v['response']]
        chosen=max(complete,key=lambda v:(v['level'],v['augmentation'])) if complete else levels[-1]
        atom={'symbol':a,'z':Z[a],'levels':levels,'selected':chosen,'basis_check':None,'diffuse_check':None}
        prev=next((v for v in levels if v['level']==chosen['level']-1 and v['augmentation']==chosen['augmentation'] and 'ccsdt' in v['response']),None)
        sparse=next((v for v in levels if v['level']==chosen['level'] and v['augmentation']==chosen['augmentation']-1 and 'ccsdt' in v['response']),None)
        for name,older in [('basis_check',prev),('diffuse_check',sparse)]:
            if older and 'ccsdt' in chosen['response']:
                atom[name]={'response_percent':100*(chosen['response']['ccsdt']['average']/older['response']['ccsdt']['average']-1),
                            'response_components_percent':{b:100*(chosen['response']['ccsdt'][b]['value']/older['response']['ccsdt'][b]['value']-1) for b in ('parallel','perpendicular')},
                            'mean_radius_percent':100*(chosen['mean_pm']/older['mean_pm']-1),
                            'removal_eV_difference':chosen['removal_eV']-older['removal_eV'] if chosen['removal_eV'] is not None else None}
        comparisons=[atom['basis_check'],atom['diffuse_check']]
        average_pass=all(c is not None and abs(c['response_percent'])<1 for c in comparisons)
        components_pass=all(c is not None and all(abs(x)<1 for x in c['response_components_percent'].values()) for c in comparisons)
        atom['response_convergence']={'threshold_percent':1,'average_pass':average_pass,
            'all_components_pass':components_pass,'passed':average_pass and components_pass,
            'scope':'Last cardinal and successive diffuse changes, each angular branch as well as the average; not a physical error bound.'}
        if a in EXPERIMENT:
            atom['benchmarks']={'removal':{'type':'NIST evaluated measured ground-level interval',
                  'value_eV':EXPERIMENT[a],'uncertainty_eV':EXPERIMENT_UNC[a],'source':NIST,
                  'residual_eV':chosen['removal_eV']-EXPERIMENT[a] if chosen['removal_eV'] else None},
                  'response':{'type':'independent nonrelativistic CCSD(T) calculation, not experiment',
                  'value':THEORY[a],'parallel':THEORY_COMPONENTS[a][0],'perpendicular':THEORY_COMPONENTS[a][1],
                  'source':PAPER,'reported_digits_only':True}}
        atoms[a]=atom
    fci_records=[]
    for p in sorted(controls.glob('*.json')):
        r=json.loads(p.read_text());fci_records.append(r);paths.append(p)
        check(r['source_sha256']==hashlib.sha256((source/'correlated_controls.py').read_bytes()).hexdigest(),'FCI solver source '+p.name)
        check(r['fci_converged'] and r['cc_converged'],'small-space FCI '+p.name)
        spin=SPINS[r['symbol']][r['charge']]/2
        check(abs(r['fci_s2']-spin*(spin+1))<1e-7,'FCI sector '+p.name)
    for a in ('C','N','O'):
        if a not in atoms:continue
        q=[next((r for r in fci_records if r['symbol']==a and r['charge']==i),None) for i in (0,1)]
        if all(q):atoms[a]['small_basis_correlation_check']={'basis':'cc-pVDZ (no augmentation)',
            'ccsdt_removal_minus_fci_eV':(q[1]['ccsdt_minus_fci']-q[0]['ccsdt_minus_fci'])*units['energy_eV'],
            'interpretation':'Finite-space correlation truncation check; not a correction to the larger basis.'}
    for a,reference in [('H',4.5),('He',original_controls['helium_infinite']['response']['polarizability'])]:
        if a in atoms and 'ccsdt' in atoms[a]['selected']['response']:
            atoms[a]['response_control']={'reference':reference,'calculated':atoms[a]['selected']['response']['ccsdt']['average'],
                'relative_percent':100*(atoms[a]['selected']['response']['ccsdt']['average']/reference-1),
                'model':'infinite-inertia nonrelativistic; H exact, He high-accuracy correlated control'}
    if strict:
        check(len(atoms)==5,'all control and target species')
        check(len(fci_records)==6,'six independent small-basis FCI checks')
        for a in ('C','N','O'):
            atom=atoms.get(a,{})
            check(atom.get('selected',{}).get('level')==4,'QZ target '+a)
            check(atom.get('selected',{}).get('augmentation',0)>=2,'diffuse target '+a)
            check(atom.get('basis_check') is not None,'basis comparison '+a)
            check(atom.get('diffuse_check') is not None,'diffuse comparison '+a)
            check(atom.get('selected',{}).get('removal_eV') is not None,'removal calculation '+a)
    archive_plan={}
    for a,atom in atoms.items():
        arrays=[p for p in statepaths if p.name.startswith(a+'-')]
        split=sum(p.stat().st_size for p in arrays)>90*1024**2
        groups=sorted({v['augmentation'] for v in atom['levels']}) if split else [None]
        atom['state_archives']=[]
        for group in groups:
            name='correlated_atoms_states-'+a+(f'-a{group}' if group is not None else '')+'.zip'
            archive_plan[name]=(a,group);atom['state_archives'].append(name)
            if group is None or group==atom['selected']['augmentation']:atom['state_archive']=name
    report={'schema':'poams.correlated_atomic_preparation.v1','model':'NR infinite-inertia point compact; ROHF-seeded UCCSD(T); all entries correlated',
            'units':units,'atoms':atoms,'validation':{'records':len(records),'state_archives':len(statepaths),
            'fci_checks':len(fci_records),'errors':errors,'passed':not errors},
            'molecules_status':'ON HOLD. Atomic controls do not certify molecular gradients, BSSE, or long optimization.',
            'density_scope':'CCSD left/right one-body density expectations; not relaxed CCSD(T) derivatives.',
            'solver_source_sha256':sorted(source_versions),
            'sources':{'NIST_ionization':NIST,'independent_response':PAPER,'cc_implementation':'https://pyscf.org/user/cc.html'}}
    out.mkdir(parents=True,exist_ok=True)
    dump(out/'correlated_atoms_results.json',report)
    if errors:raise RuntimeError('\n'.join(errors))
    (out/'correlated_atoms_data.js').write_text('window.POAMS_CORRELATED_ATOMS = '+json.dumps(report,separators=(',',':'),allow_nan=False)+';\n')
    with zipfile.ZipFile(out/'correlated_atoms_records.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:archive.write(path,'records/'+path.name)
        for name in ['correlated_atoms.py','correlated_campaign.py','correlated_controls.py','correlated_atoms_build.py',
                     'correlated_atoms_notes.md','whole_atom_handoff.json']:
            archive.write(source/name,'scripts/'+name)
        for path in sorted((source/'correlated_atoms_sources').glob('*.py')):
            archive.write(path,'scripts/correlated_atoms_sources/'+path.name)
        archive.writestr('HANDOFF.json',json.dumps(report,indent=2))
    for name,(a,group) in archive_plan.items():
        with zipfile.ZipFile(out/name,'w',compression=zipfile.ZIP_DEFLATED) as archive:
            for path in statepaths:
                if path.name.startswith(a+'-') and (group is None or f'-a{group}-' in path.name):archive.write(path,'states/'+path.name)
            for path in paths:
                if path.name.startswith(a+'-') and (group is None or f'-a{group}-' in path.name):archive.write(path,'records/'+path.name)
        if (out/name).stat().st_size>90*1024**2:raise RuntimeError('State archive requires further splitting: '+name)
    print(json.dumps(report['validation']))
    for a,r in atoms.items():
        s=r['selected'];print(a,s['basis'],s['removal_eV'],s['response'].get('ccsdt',{}).get('average'),r['basis_check'],r['diffuse_check'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--roots',nargs='+',type=Path,required=True)
    p.add_argument('--controls',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--strict',action='store_true');a=p.parse_args()
    build(a.roots,a.controls,a.out,a.strict)
