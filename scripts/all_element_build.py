#!/usr/bin/env python3
"""Build auditable display data from saved atomic calculations; no new solve."""
import argparse
import csv
import io
import hashlib
import json
import zipfile
from pathlib import Path
import numpy as np

def state_archives(states,rows,dest,max_bytes=40*1024*1024):
    """Keep reusable state downloads below ordinary repository file limits."""
    index_path=dest/'all_element_state_index.json'
    previous=json.loads(index_path.read_text())['archives'] if index_path.exists() else []
    chunks=[];chunk=[];size=0
    for file,row in zip(states,rows):
        if not file.is_file():raise FileNotFoundError(file)
        length=file.stat().st_size
        if length>max_bytes:raise ValueError(f'State exceeds archive budget: {file}')
        if chunk and size+length>max_bytes:
            chunks.append(chunk);chunk=[];size=0
        chunk.append((file,row));size+=length
    if chunk:chunks.append(chunk)
    archives=[]
    for chunk in chunks:
        name=f"all_element_states_{chunk[0][1]['z']:03d}-{chunk[-1][1]['z']:03d}.zip"
        path=dest/name;temporary=path.with_suffix('.tmp')
        members=[]
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_STORED) as archive:
            for file,row in chunk:
                archive.write(file,file.name)
                row['state_archive']=name;row['state_file']=file.name
                members.append({'z':row['z'],'symbol':row['symbol'],'level':row['level'],
                                'file':file.name,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
        temporary.replace(path)
        archives.append({'file':name,'bytes':path.stat().st_size,
                         'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'members':members})
    index={'schema':'poams.atomic.state-archives.v1','states':len(rows),
           'scope':'Selected highest validated basis per element; NumPy allow_pickle=False. See all_element_notes.md.',
           'archives':archives}
    index_path.write_text(json.dumps(index,indent=2)+'\n')
    current={a['file'] for a in archives}
    for old in previous:
        name=old['file']
        if name not in current and Path(name).name==name and name.startswith('all_element_states_') and name.endswith('.zip'):
            (dest/name).unlink(missing_ok=True)
    return archives


def build(source,dest):
    source=Path(source);dest=Path(dest)
    units=json.loads((dest/'whole_atom_handoff.json').read_text())['units']
    # Existing atomic units are the same measured-alpha calibration.
    a=5.291772105467402e-11;energy=4.359744722158896e-18;rate=6.579683920428908e15
    data={'schema':'poams.atomic.meanfield.v1','date':'2026-09-15',
          'method':'Configuration-seeded UHF; scalar X2C; all-electron Dyall + diffuse functions',
          'units':{'length_m':a,'energy_J':energy,'rate_Hz':rate,'energy_eV':27.21138624568717},
          'scope':'Neutral, fixed point compact component; external reference only. Not correlated precision or isotope-core structure.',
          'elements':{},'counts':{}}
    diagnostic_path=dest/'all_element_response_diagnostics.json'
    diagnostics=json.loads(diagnostic_path.read_text())['states'] if diagnostic_path.exists() else []
    full=[];states=[];failures=[]
    for z in range(1,119):
        levels={}
        paths={}
        for level in [2,3]:
            files=list(source.glob(f'{z:03d}-*-v{level}z.json'))
            if files:
                raw=json.loads(files[0].read_text());levels[level]=raw;paths[level]=files[0]
        valid=[level for level,r in levels.items() if r['validated']]
        if not valid:
            failures.append(z);continue
        level=max(valid);r=levels[level]
        state_path=source/(paths[level].stem+'-state.npz')
        with np.load(state_path,allow_pickle=False) as saved:
            nmo=int(saved['mo_coeff'].shape[2])
        response=r['response']; chi=response.get('relaxed_isotropic') if response.get('certified') else None
        comparison={}
        if 2 in levels and 3 in levels and levels[2]['validated'] and levels[3]['validated']:
            low=levels[2];high=levels[3]
            comparison={'energy_difference':high['energy']-low['energy'],
                        'mean_radius_percent':100*(high['radial']['mean_r']/low['radial']['mean_r']-1),
                        'response_percent':None}
            if low['response'].get('certified') and high['response'].get('certified'):
                comparison['response_percent']=100*(high['response']['relaxed_isotropic']/low['response']['relaxed_isotropic']-1)
        angular=[sum(o['l']==l for o in r['occupied_orbitals']) for l in range(4)]
        row={k:r[k] for k in ['z','symbol','level','energy','one_body_energy','direct_energy','exchange_energy',
                              'gradient_norm','integral_energy_difference','direct_gradient_norm','input_spin_2Ms',
                              'spin_squared','orbital_L_squared','koopmans_removal','nao','occupied_orbitals']}
        row.update({'nmo_per_spin':nmo,'overlap_eigenvalue_cutoff':1e-6,
                    'calculated_angular_counts':angular,'input_angular_counts':r['input_configuration_l_counts'],
                    'energy_rate_PHz':r['energy']*rate/1e15,
                    'mean_pm':r['radial']['mean_r']*a*1e12,'rms_pm':r['radial']['rms_r']*a*1e12,
                    'r90_pm':r['radial']['r90']*a*1e12,
                    'radial_profile':[[round(x*a*1e12,9),float(f'{p/(a*1e12):.10g}')] for x,p in r['radial']['profile']],
                    'radial_normalization_error':r['radial']['normalization']-z,
                    'response':{'status':response['status'],'relaxed':chi,
                                'tensor':response.get('relaxed_tensor') if chi is not None else None,
                                'frozen':response.get('frozen_isotropic'),
                                'frozen_C6':response.get('frozen_C6_self'),
                                'scaled_residual':response.get('scaled_linear_residual')},
                    'comparison':comparison,'levels':{},
                    'state_scope':'Stationary determinant in the specified spin sector, not a global ground-term proof.'})
        # Attach a diagnosis only to the exact saved state it tested.
        for diagnostic in diagnostics:
            if (diagnostic['z']==z and diagnostic['level']==level and
                diagnostic['source_sha256']==hashlib.sha256(paths[level].read_bytes()).hexdigest() and
                diagnostic['source_state_sha256']==hashlib.sha256(state_path.read_bytes()).hexdigest()):
                row['response']['stability_diagnosis']='Energy-lowering orbital rotation confirmed: this stationary reference is unstable in a tested dipole-coupled mode.'
        ion=paths[level].with_name(paths[level].stem+'-removal.json')
        row['removal']=json.loads(ion.read_text()) if ion.exists() else None
        for lv,rr in levels.items():
            row['levels'][str(lv)]={'validated':rr['validated'],'energy':rr['energy'],
                                  'mean_pm':rr['radial']['mean_r']*a*1e12,
                                  'response':rr['response'].get('relaxed_isotropic') if rr['response'].get('certified') else None}
        data['elements'][str(z)]=row
        full.append((paths[level].name,r));states.append(state_path)
    rows=list(data['elements'].values())
    data['counts']={'states':len(rows),'two_basis':sum(bool(r['comparison']) for r in rows),
                    'relaxed_responses':sum(r['response']['relaxed'] is not None for r in rows),
                    'relaxed_removals':sum(bool(r['removal'] and r['removal'].get('certified')) for r in rows),
                    'basis_sensitive_response':sum(abs(r['comparison'].get('response_percent') or 0)>5 for r in rows),
                    'direct_gradient_flagged':sum(r['direct_gradient_norm']>=3e-5 for r in rows),
                    'missing_states':failures}
    data['state_archives']=state_archives(states,rows,dest)
    (dest/'all_element_results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    (dest/'all_element_data.js').write_text('// Generated by all_element_build.py.\nwindow.POAMS_ALL_ELEMENTS = '+json.dumps(data,separators=(',',':'),allow_nan=False)+';\n')
    stream=io.StringIO();writer=csv.writer(stream,lineterminator='\n')
    writer.writerow(['Z','symbol','basis_level','external_E_over_h_PHz','mean_envelope_pm','rms_envelope_pm','r90_pm','relaxed_removal_PHz',
                     'relaxed_static_chi','frozen_static_chi','basis_chi_difference_percent','spin_squared','orbital_L_squared',
                     'specified_2Ms','SCF_gradient','direct_integral_energy_difference_Estar'])
    for r in rows:
        writer.writerow([r['z'],r['symbol'],r['level'],r['energy_rate_PHz'],r['mean_pm'],r['rms_pm'],r['r90_pm'],
                         r['removal']['interval']*rate/1e15 if r['removal'] and r['removal'].get('certified') else None,
                         r['response']['relaxed'],r['response']['frozen'],r['comparison'].get('response_percent'),
                         r['spin_squared'],r['orbital_L_squared'],r['input_spin_2Ms'],r['gradient_norm'],r['integral_energy_difference']])
    (dest/'all_element_table.csv').write_text(stream.getvalue())
    with zipfile.ZipFile(dest/'all_element_records.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for name,r in full:archive.write(source/name,name)
        for f in source.glob('*-v?z.json'):
            if f.name not in {name for name,_ in full}:archive.write(f,'comparison/'+f.name)
        for f in source.glob('*-removal.json'):archive.write(f,'removal/'+f.name)
        for name in ['all_element_atoms.py','all_element_campaign.py','all_element_build.py','all_element_verify.py','all_element_notes.md','all_element_response_diagnostics.py','all_element_response_diagnostics.json']:
            if (dest/name).exists():archive.write(dest/name,name)
    print(json.dumps(data['counts']))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True)
    p.add_argument('--dest',default=str(Path(__file__).resolve().parent));a=p.parse_args();build(a.source,a.dest)
