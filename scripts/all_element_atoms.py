#!/usr/bin/env python3
"""Neutral-atom mean-field campaign, Z=1..118 (not a correlated-core solver).

Requires Python, numpy, scipy, PySCF 2.14.0. See all_element_notes.md.
Example: python all_element_atoms.py --z 6 --level 2 --out RESULTS
Units inside the solver: calibrated a_*=hbar/(m*c*alpha), E_*=m*c^2*alpha^2.
The scalar X2C/UHF effective operator is conventional-equivalent mathematics;
neither it nor the imported Gaussian representation is derived from ontology.
"""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import os
import time
from pathlib import Path

os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.integrate import simpson, cumulative_simpson
from scipy.sparse.linalg import LinearOperator, gmres, minres
from pyscf import gto, scf, df
from pyscf.data import elements
from pyscf.lib import param
from pyscf.scf import ucphf, atom_hf, _response_functions  # registers gen_response
import pyscf

ALPHA = 0.0072973525643
param.LIGHT_SPEED = 1 / ALPHA
C = param.LIGHT_SPEED


def spin_sector(z):
    """Input sector: aligned open-shell spins of the tabulated configuration."""
    return int(sum(min(n % (2*(2*l+1)), 2*(2*l+1)-n % (2*(2*l+1)))
                   for l, n in enumerate(elements.CONFIGURATION[z])))


def basis_for(z, level):
    """Dyall vNz plus two geometrically continued diffuse functions per l.

    Diffuse ratio 1/3 is fixed before results; no response target is fitted.
    Include l through the highest occupied l + 1 (at least p).
    """
    sym = elements.ELEMENTS[z]
    basis = copy.deepcopy(gto.basis.load(f'dyall-v{level}z', sym))
    max_occ = max(l for l, n in enumerate(elements.CONFIGURATION[z]) if n)
    for l in range(max(1, max_occ+1)+1):
        exps = [float(row[0]) for shell in basis if shell[0] == l
                for row in shell[2:]]
        if not exps:
            # A missing response channel gets the smallest occupied-l exponent.
            exps = [float(row[0]) for shell in basis if shell[0] == max_occ
                    for row in shell[2:]]
        low = min(exps)
        basis += [[l, 0, [low/3, 1.0]], [l, 0, [low/9, 1.0]]]
    return basis


def make_mf(z, level, verbose=0, direct=False, charge=0, spin=None):
    mol = gto.M(atom=f'{elements.ELEMENTS[z]} 0 0 0', unit='Bohr',
                basis=basis_for(z, level), spin=spin_sector(z) if spin is None else spin,
                charge=charge, verbose=verbose, max_memory=2200)
    mf = scf.UHF(mol).sfx2c1e()
    if not direct:
        mf = mf.density_fit(auxbasis=df.addons.aug_etb(mol,beta=2.0))
    mf.init_guess = '1e'
    # Heavy-atom absolute energies cannot use a sub-roundoff stopping test.
    # The independent orbital-gradient criterion remains unchanged.
    mf.conv_tol = max(2e-10,1e-12*z*z)
    mf.conv_tol_grad = 2e-6
    mf.max_cycle = 180
    mf.chkfile = None
    return mf


def configuration_seed(mf,z):
    """Spherical radial optimization, then integer Hund slots as a UHF seed.

    Hydrogenic core guesses can trap neutral Li in 1s2 2p, for example.
    The supplied angular configuration is an input, not a derived prediction.
    Subsequent UHF relaxation is not constrained to these orbital shapes.
    """
    if z==1:
        return None,{'method':'one-entry direct diagonalization'}
    seed=atom_hf.AtomSphAverageRHF(mf.mol).sfx2c1e()
    if getattr(mf,'with_df',None) is not None:
        seed=seed.density_fit(auxbasis=mf.with_df.auxbasis)
    seed.atomic_configuration=elements.CONFIGURATION
    seed.init_guess='1e';seed.conv_tol=1e-8;seed.max_cycle=140
    seed.kernel()
    c=seed.mo_coeff
    occupation=np.zeros((2,c.shape[1]));offset=0
    mol=mf.mol
    for l in sorted(set(angular_labels(mol))):
        g=2*l+1
        nr=sum(mol.bas_nctr(i) for i in range(mol.nbas) if mol.bas_angular(i)==l)
        total=elements.CONFIGURATION[z][l] if l<4 else 0
        complete,partial=divmod(total,2*g)
        occupation[:,offset:offset+complete*g]=1
        start=offset+complete*g
        occupation[0,start:start+min(g,partial)]=1
        occupation[1,start:start+max(0,partial-g)]=1
        offset+=nr*g
    assert np.sum(occupation)==z
    assert np.sum(occupation[0]-occupation[1])==mf.mol.spin
    return np.array([(c*occ)@c.T for occ in occupation]),{
        'method':'configuration-seeded spherical radial reference / integer Hund slots',
        'radial_seed_converged':bool(seed.converged),'radial_seed_energy':float(seed.e_tot)}


def solve(mf, dm0=None):
    mf.kernel(dm0=dm0)
    attempts = [{'method': 'DIIS', 'converged': bool(mf.converged),
                 'cycles': int(mf.cycles), 'energy': float(mf.e_tot)}]
    if not mf.converged:
        mf.diis_space = 16
        mf.level_shift = .2
        mf.max_cycle = 160
        mf.kernel(dm0=mf.make_rdm1())
        attempts.append({'method': 'shifted DIIS', 'converged': bool(mf.converged),
                         'cycles': int(mf.cycles), 'energy': float(mf.e_tot)})
        mf.level_shift = 0
        if mf.converged:
            mf.kernel(dm0=mf.make_rdm1())
    if not mf.converged:
        new = mf.newton()
        new.max_cycle = 80
        new.kernel(mf.mo_coeff, mf.mo_occ)
        attempts.append({'method': 'second order', 'converged': bool(new.converged),
                         'energy': float(new.e_tot)})
        mf.mo_coeff, mf.mo_energy, mf.mo_occ = new.mo_coeff, new.mo_energy, new.mo_occ
        mf.e_tot, mf.converged = new.e_tot, new.converged
    return attempts


def angular_labels(mol):
    ls = np.zeros(mol.nao, dtype=int)
    loc = mol.ao_loc_nr()
    for i in range(mol.nbas):
        ls[loc[i]:loc[i+1]] = mol.bas_angular(i)
    return ls


def spherical_average(d, mol):
    """Exact rotation average of a one-centre real-spherical AO density."""
    ls = angular_labels(mol)
    out = np.zeros_like(d)
    for l in np.unique(ls):
        idx = np.flatnonzero(ls == l)
        g = 2*l+1
        n = len(idx)//g
        block = d[np.ix_(idx, idx)].reshape(n,g,n,g)
        avg = np.einsum('imjm->ij', block)/g
        out[np.ix_(idx, idx)] = np.einsum('ij,mn->imjn', avg, np.eye(g)).reshape(n*g,n*g)
    return out


def radial_density(mf, samples=2401):
    """Physical large+small component radial density in spin-free X2C.

    Small components use restricted kinetic balance: sigma.p/(2c) X R C.
    The scalar spin trace leaves the sum of squared spatial gradients.
    Rotation averaging commutes with the scalar atomic X2C transform.
    """
    mol = mf.mol
    d = mf.make_rdm1().sum(axis=0)
    x = mf.with_x2c.get_xmat()
    rmat = mf.with_x2c._get_rmat(x)
    large = rmat @ d @ rmat.T
    small = x @ large @ x.T
    ov = mol.intor_symmetric('int1e_ovlp')
    kin = mol.intor_symmetric('int1e_kin')
    norm_matrix = np.einsum('ij,ji', large, ov) + np.einsum('ij,ji', small, kin)/(2*C*C)
    large = spherical_average(large, mol)
    small = spherical_average(small, mol)
    radius = np.geomspace(1e-9, 300., samples)
    xyz = np.zeros((samples,3)); xyz[:,2] = radius
    ao = mol.eval_gto('GTOval_sph_deriv1', xyz)
    rho = np.einsum('pi,ij,pj->p', ao[0],large,ao[0],optimize=True)
    rho += sum(np.einsum('pi,ij,pj->p',a,small,a,optimize=True)
               for a in ao[1:])/(4*C*C)
    probability = 4*np.pi*radius**2*rho
    logr = np.log(radius)
    moments = [float(simpson(probability*radius**(k+1), x=logr)) for k in range(3)]
    cumulative = cumulative_simpson(probability*radius, x=logr, initial=0.)
    r90 = float(np.interp(.9*mol.nelectron, cumulative, radius))
    # An independent nested-grid check, not an arbitrary normalization repair.
    halfnorm = float(simpson((probability*radius)[::2],x=logr[::2]))
    halfmean = float(simpson((probability*radius**2)[::2],x=logr[::2]))/mol.nelectron
    keep = np.unique(np.r_[np.arange(0,samples,6), samples-1])
    return {'normalization': moments[0], 'matrix_normalization':float(norm_matrix),
            'coarse_normalization':halfnorm,
            'mean_r':moments[1]/mol.nelectron,
            'rms_r':float(np.sqrt(moments[2]/mol.nelectron)), 'r90':r90,
            'mean_grid_difference':abs(halfmean-moments[1]/mol.nelectron),
            'profile':np.column_stack((radius[keep],probability[keep]/mol.nelectron)).tolist()}


def dipole_matrix(mf):
    mol = mf.mol
    prp = mol.intor_symmetric('int1e_sprsp').reshape(3,4,mol.nao,mol.nao)[:,3]
    return mf.with_x2c.picture_change(('int1e_r', prp/(4*C*C)))


def response(mf, relaxed=True):
    """Static tensor in the fixed projected-X2C model; no force taxonomy."""
    dip = dipole_matrix(mf)
    blocks=[]; gaps=[]; cv=[]; ci=[]; parity_masks=[]
    ov=mf.get_ovlp();ls=angular_labels(mf.mol)
    parity_ao=ov*((-1.)**ls)[None,:]
    for e,c,o in zip(mf.mo_energy,mf.mo_coeff,mf.mo_occ):
        occ=o>0; cv.append(c[:,~occ]); ci.append(c[:,occ])
        blocks.append(np.einsum('pa,xpq,qi->xai',cv[-1],dip,ci[-1],optimize=True))
        gaps.append(e[~occ,None]-e[occ])
        parity=np.einsum('pi,pq,qi->i',c,parity_ao,c,optimize=True)
        mask=(parity[~occ,None]*parity[occ]<0) if np.min(np.abs(parity))>.99999 else np.ones_like(gaps[-1],dtype=bool)
        parity_masks.append(mask)
        blocks[-1]*=mask
    bad = [(gap < 1e-7)&(np.max(np.abs(b),axis=0)>1e-7) for gap,b in zip(gaps,blocks)]
    if any(np.any(b) for b in bad):
        return {'status':'withheld: nonpositive or unresolved dipole-coupled gap'}
    safe=[np.where(g>1e-7,g,np.inf) for g in gaps]
    frozen=sum(2*np.einsum('xai,yai,ai->xy',b,b,1/g,optimize=True)
               for b,g in zip(blocks,safe))
    vals=[]; strengths=[]
    for g,b in zip(gaps,blocks):
        ss=np.sum(b*b,axis=0)/3
        select=(g>1e-7)&(ss>1e-13)
        vals.extend(g[select]); strengths.extend(ss[select])
    vals=np.array(vals); strengths=np.array(strengths)
    # Independent spectral and tensor definitions must agree.
    spectral_static=float(2*np.sum(strengths/vals))
    xx,ww=np.polynomial.legendre.leggauss(96)
    theta=(xx+1)*np.pi/4
    omega=np.tan(theta); weights=ww*np.pi/4/np.cos(theta)**2
    imaginary=2*np.sum(strengths[:,None]*vals[:,None]/(vals[:,None]**2+omega**2),axis=0)
    c6=float(3/np.pi*np.dot(weights,imaginary**2))
    result={'status':'frozen-reference only', 'frozen_tensor':frozen.tolist(),
            'frozen_isotropic':float(np.trace(frozen)/3),
            'frozen_spectral_static':spectral_static,'frozen_C6_self':c6,
            'transition_energies':vals.tolist(),'transition_strengths':strengths.tolist(),
            'minimum_dipole_gap':float(min(vals))}
    if mf.mol.nelectron == 1:
        # PySCF's one-entry solver diagonalizes h, not an N-entry Fock map.
        # There is no pair response kernel: adding one would create self-action.
        result.update({'status':'one-entry finite-basis response',
                       'relaxed_tensor':frozen.tolist(),
                       'relaxed_isotropic':float(np.trace(frozen)/3),
                       'linear_residual':0.,'scaled_linear_residual':0.,
                       'tensor_asymmetry':0.,'certified':True})
        return result
    if not relaxed:
        return result
    shapes=[b.shape[1:] for b in blocks]
    split=int(np.prod(shapes[0]))
    vresp=mf.gen_response(hermi=1)
    def vind(flat):
        nset=len(flat)
        us=[flat[:,:split].reshape((nset,)+shapes[0]),flat[:,split:].reshape((nset,)+shapes[1])]
        ds=[]
        for u,v,i,mask in zip(us,cv,ci,parity_masks):
            q=np.einsum('pa,xai,qi->xpq',v,u*mask,i,optimize=True)
            ds.append(q+q.transpose(0,2,1))
        veff=vresp(np.array(ds))
        return np.hstack([(np.einsum('pa,xpq,qi->xai',v,vv,i,optimize=True)*mask).reshape(nset,-1)
                          for v,vv,i,mask in zip(cv,veff,ci,parity_masks)])
    try:
        rhs=np.hstack([b.reshape(3,-1) for b in blocks])
        gapflat=np.hstack([g.ravel() for g in gaps])
        nvar=len(gapflat)
        operator=LinearOperator((nvar,nvar),matvec=lambda u:gapflat*u+vind(u[None])[0])
        preconditioner=LinearOperator((nvar,nvar),matvec=lambda u:u/(np.abs(gapflat)+.2))
        # MINRES exploits the symmetric response Hessian, including indefinite
        # trial states. The shift is only in the positive preconditioner.
        # Acceptance is still based on an independently recomputed residual.
        flat=np.zeros_like(rhs);minres_info=[];iteration_counts=[]
        for axis in range(3):
            counter=[0]
            def count_iteration(vector):counter[0]+=1
            flat[axis],info=minres(operator,-rhs[axis],M=preconditioner,
                                  rtol=1e-11,maxiter=300,callback=count_iteration)
            minres_info.append(int(info));iteration_counts.append(counter[0])
        residual_vector=flat*gapflat+vind(flat)+rhs
        gmres_info=[]
        if np.max(np.abs(residual_vector)/(np.abs(gapflat)+.2)) > 2e-9:
            for axis in range(3):
                flat[axis],info=gmres(operator,-rhs[axis],x0=flat[axis],M=preconditioner,
                                      rtol=1e-10,atol=1e-11,restart=35,maxiter=8)
                gmres_info.append(int(info))
            residual_vector=flat*gapflat+vind(flat)+rhs
        us=[flat[:,:split].reshape((3,)+shapes[0]),flat[:,split:].reshape((3,)+shapes[1])]
        residual=float(np.max(np.abs(residual_vector)))
        scaled_residual=float(np.max(np.abs(residual_vector)/(np.abs(gapflat)+.2)))
        tensor=-2*sum(np.einsum('xai,yai->xy',b,u) for b,u in zip(blocks,us))
        asym=float(np.max(np.abs(tensor-tensor.T)))
        tensor=(tensor+tensor.T)/2
        stable=bool(np.min(np.linalg.eigvalsh(tensor))>0 and scaled_residual<2e-8 and asym<2e-5)
        result.update({'status':'relaxed finite-basis response' if stable else 'relaxed response not certified',
                       'relaxed_tensor':tensor.tolist(), 'relaxed_isotropic':float(np.trace(tensor)/3),
                       'linear_residual':residual,'scaled_linear_residual':scaled_residual,
                       'tensor_asymmetry':asym,'minres_info':minres_info,
                       'response_iterations':iteration_counts,'gmres_info':gmres_info,'certified':stable})
    except Exception as exc:
        result['relaxed_failure']=str(exc)
    return result


def calculate(z,level,out,do_response=True,direct=False):
    folder=Path(out);folder.mkdir(parents=True,exist_ok=True)
    with (folder/f'{z:03d}-v{level}z.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _calculate_locked(z,level,out,do_response,direct)


def _calculate_locked(z,level,out,do_response=True,direct=False):
    started=time.time()
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    name=f'{z:03d}-{elements.ELEMENTS[z]}-v{level}z'
    final=out/f'{name}.json'
    if final.exists():
        print('EXISTS',name,flush=True); return
    progress=out/f'{name}-progress.json'
    def stage(label):
        progress.write_text(json.dumps({'z':z,'level':level,'stage':label,'started':started,'updated':time.time()}))
    checkpoint=out/f'{name}-checkpoint.json'
    statefile=out/f'{name}-state.npz'
    if checkpoint.exists() and statefile.exists():
        mf,result=load_reference(checkpoint,statefile)
        return finish_reference(mf,result,final,progress,started,checkpoint,do_response)
    stage('self-consistent state')
    mf=make_mf(z,level,direct=direct)
    dm0=None
    seed=out/f'{z:03d}-{elements.ELEMENTS[z]}-v2z-state.npz'
    if level==3 and seed.exists():
        with np.load(seed,allow_pickle=False) as saved:
            old=gto.loads(str(saved['mol_json']))
            cc=saved['mo_coeff'];oo=saved['mo_occ']
            olddm=np.array([(c*o)@c.T for c,o in zip(cc,oo)])
            dm0=scf.addons.project_dm_nr2nr(old,olddm,mf.mol)
    if dm0 is None:
        dm0,seed_info=configuration_seed(mf,z)
    else:
        seed_info={'method':'projected validated double-zeta state'}
    attempts=solve(mf,dm0)
    stage('state diagnostics')
    m=mf.mol;dm=mf.make_rdm1();d=dm.sum(axis=0)
    ov=mf.get_ovlp();h=mf.get_hcore();f=mf.get_fock(dm=dm)
    grad=float(np.linalg.norm(mf.get_grad(mf.mo_coeff,mf.mo_occ,f)))
    orth=max(float(np.max(np.abs(c.T@ov@c-np.eye(c.shape[1])))) for c in mf.mo_coeff)
    count=float(np.einsum('ij,ji',d,ov))
    one=float(np.einsum('ij,ji',d,h))
    j,k=mf.get_jk(m,dm)
    direct_term=float(.5*np.einsum('ij,ji',d,j.sum(axis=0)))
    exchange=float(-.5*np.einsum('sij,sji',dm,k))
    # Independent direct-integral check at this density; not a fitted datum.
    if not direct:
        raw=scf.UHF(m)
        jd,kd=raw.get_jk(m,dm)
        direct_integral_energy=one+.5*np.einsum('ij,ji',d,jd.sum(axis=0))-.5*np.einsum('sij,sji',dm,kd)
        direct_gradient=float(np.linalg.norm(mf.get_grad(mf.mo_coeff,mf.mo_occ,h+jd.sum(axis=0)-kd)))
    else:
        direct_integral_energy=one+direct_term+exchange
        direct_gradient=grad
    ls=angular_labels(m)
    lang=-1j*m.intor('int1e_cg_irxp',comp=3)
    l2ao=ov*(ls*(ls+1))[None,:]
    lmean=np.einsum('ij,xji->x',d,lang).real
    # Explicit matrix products avoid numpy's six-index uncontracted path.
    exchange_l2=0j
    for ds in dm:
        for generator in lang:
            q=ds@generator
            exchange_l2+=np.einsum('ij,ji',q,q)
    l2=(np.einsum('ij,ji',d,l2ao)-exchange_l2+np.dot(lmean,lmean)).real
    orbital=[]
    for s,(en,co,oc) in enumerate(zip(mf.mo_energy,mf.mo_coeff,mf.mo_occ)):
        for idx in np.flatnonzero(oc>0):
            c=co[:,idx];weights=[]
            for l in np.unique(ls):
                ix=np.flatnonzero(ls==l)
                weights.append((int(l),float(c[ix]@ov[np.ix_(ix,ix)]@c[ix])))
            l,w=max(weights,key=lambda q:q[1])
            orbital.append({'spin':s,'index':int(idx),'energy':float(en[idx]),
                            'l':l,'l_weight':w})
    result={'z':z,'symbol':elements.ELEMENTS[z],'level':level,'pyscf':pyscf.__version__,
            'basis':f'dyall-v{level}z + 2 diffuse/channel', 'alpha':ALPHA,
            'integrals':'direct' if direct is True else 'RI-JK, generated even-tempered auxiliary beta=2',
            'model':'UHF / spin-free one-electron X2C / point fixed compact component',
            'input_spin_2Ms':spin_sector(z),'input_configuration_l_counts':elements.CONFIGURATION[z],
            'nao':m.nao, 'energy':float(mf.e_tot),'one_body_energy':one,
            'direct_energy':direct_term,'exchange_energy':exchange,
            'energy_sum_error':abs(one+direct_term+exchange-mf.e_tot),
            'direct_integral_energy_at_state':float(direct_integral_energy),
            'integral_energy_difference':float(direct_integral_energy-mf.e_tot),
            'direct_gradient_norm':direct_gradient,
            'converged':bool(mf.converged),'gradient_norm':grad,'count':count,
            'orthonormality_error':orth,'attempts':attempts,'seed':seed_info,
            'spin_squared':float(mf.spin_square()[0]),'orbital_L_squared':float(l2),
            'occupied_orbitals':orbital,
            'koopmans_removal':float(-max(row['energy'] for row in orbital)),
            'radial':radial_density(mf)}
    result['validated']=bool(mf.converged and grad<3e-5 and orth<2e-7
                             and abs(count-z)<2e-7
                             and abs(result['radial']['normalization']-z)<2e-5
                             and result['energy_sum_error']<1e-7)
    # Save the costly neutral state before entering a potentially long response.
    # A restart can resume it without treating a pending response as complete.
    np.savez_compressed(statefile,mo_coeff=mf.mo_coeff,
                        mo_energy=mf.mo_energy,mo_occ=mf.mo_occ,
                        basis_json=json.dumps(m.basis),mol_json=m.dumps())
    result['neutral_seconds']=time.time()-started
    temp=checkpoint.with_suffix('.tmp')
    temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temp.replace(checkpoint)
    return finish_reference(mf,result,final,progress,started,checkpoint,do_response)


def finish_reference(mf,result,final,progress,started,checkpoint,do_response=True):
    progress.write_text(json.dumps({'z':result['z'],'level':result['level'],
        'stage':'directional response','started':started,'updated':time.time()}))
    if result['converged'] and result['gradient_norm']<3e-5:
        result['response']=response(mf,do_response)
    else:
        result['response']={'status':'withheld: reference not stationary'}
    result['seconds']=time.time()-started
    tmp=final.with_suffix('.tmp')
    tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    tmp.replace(final)
    progress.write_text(json.dumps({'z':result['z'],'level':result['level'],
        'stage':'finished','started':started,'updated':time.time()}))
    print(json.dumps({'z':result['z'],'symbol':result['symbol'],'level':result['level'],
                      'validated':result['validated'],'energy':result['energy'],
                      'gradient':result['gradient_norm'],'response':result['response'].get('status'),
                      'seconds':round(result['seconds'],2)}),flush=True)


def load_reference(path,statefile=None):
    path=Path(path);raw=json.loads(path.read_text())
    mf=make_mf(raw['z'],raw['level'],direct=raw.get('integrals')=='direct')
    statefile=statefile or path.with_name(path.stem+'-state.npz')
    with np.load(statefile,allow_pickle=False) as saved:
        mf.mo_coeff=saved['mo_coeff'];mf.mo_energy=saved['mo_energy'];mf.mo_occ=saved['mo_occ']
    mf.e_tot=raw['energy'];mf.converged=raw['converged']
    return mf,raw


def removal_calculation(path):
    """Separately relaxed +1 state, trying removal from each occupied spin end."""
    path=Path(path);output=path.with_name(path.stem+'-removal.json')
    if output.exists():return json.loads(output.read_text())
    started=time.time();mf,raw=load_reference(path)
    if not raw['validated']:return {'certified':False,'reason':'neutral state not validated'}
    z=raw['z'];attempts=[]
    if z==1:
        result={'z':z,'level':raw['level'],'certified':True,'interval':-raw['energy'],
                'ion_energy':0.,'ion_spin_2Ms':0,'attempts':[],'seconds':time.time()-started}
    else:
        h0=mf.get_hcore();dm=mf.make_rdm1()
        for removed_spin in [0,1]:
            if not np.any(mf.mo_occ[removed_spin]>0):continue
            ns=[int(np.sum(x)) for x in mf.mo_occ];ns[removed_spin]-=1
            spin=abs(ns[0]-ns[1]);ion=make_mf(z,raw['level'],charge=1,spin=spin)
            ion.with_df=mf.with_df
            ion.get_hcore=lambda mol=None:h0
            seed=dm.copy();indices=np.flatnonzero(mf.mo_occ[removed_spin]>0)
            j=indices[np.argmax(mf.mo_energy[removed_spin][indices])]
            c=mf.mo_coeff[removed_spin][:,j];seed[removed_spin]-=np.outer(c,c)
            if ns[0]<ns[1]:seed=seed[::-1]
            history=solve(ion,seed)
            gd=float(np.linalg.norm(ion.get_grad(ion.mo_coeff,ion.mo_occ,ion.get_fock(dm=ion.make_rdm1()))))
            attempts.append({'removed_spin':removed_spin,'ion_spin_2Ms':spin,'energy':float(ion.e_tot),
                             'converged':bool(ion.converged),'gradient_norm':gd,'history':history,
                             'spin_squared':float(ion.spin_square()[0])})
        valid=[r for r in attempts if r['converged'] and r['gradient_norm']<3e-5]
        best=min(valid,key=lambda a:a['energy']) if valid else None
        result={'z':z,'level':raw['level'],'certified':bool(best),'attempts':attempts,
                'seconds':time.time()-started,
                'scope':'Lowest converged +1 state among the two tested spin-removal seeds, not global multiplet proof.'}
        if best:result.update({'interval':best['energy']-raw['energy'],'ion_energy':best['energy'],
                               'ion_spin_2Ms':best['ion_spin_2Ms']})
    temp=output.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temp.replace(output)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--z',type=int,required=True)
    p.add_argument('--level',type=int,choices=[2,3],default=2)
    p.add_argument('--out',required=True)
    p.add_argument('--frozen-only',action='store_true')
    p.add_argument('--direct',action='store_true')
    a=p.parse_args()
    if not 1<=a.z<=118: p.error('z must be 1..118')
    calculate(a.z,a.level,a.out,not a.frozen_only,a.direct)
