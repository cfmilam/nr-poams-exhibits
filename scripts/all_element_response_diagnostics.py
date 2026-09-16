"""Saved-reference response-direction curvature audit. Does not modify campaign states."""
import os
os.environ['VECLIB_MAXIMUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,time,hashlib,argparse
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import LinearOperator,minres,gmres
sys.path.insert(0,str(Path(__file__).resolve().parent))
import all_element_atoms as a
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
parser.add_argument('--states',nargs='+',default=['39:3','40:3'],metavar='Z:LEVEL',
                    help='Saved negative-response references to diagnose (default: 39:3 40:3).')
parser.add_argument('--stable-energy-differences',action='store_true',
                    help='Also evaluate exact quadratic-functional energy increments without subtracting total energies.')
args=parser.parse_args()
source=args.source;out=args.out
try:
    targets=[tuple(map(int,token.split(':'))) for token in args.states]
    if any(len(pair)!=2 or not 1<=pair[0]<=118 or pair[1] not in (2,3) for pair in targets):
        raise ValueError
except ValueError:
    parser.error('--states must contain Z:LEVEL pairs with Z=1..118 and LEVEL=2 or 3')
out.parent.mkdir(parents=True,exist_ok=True)
result={'scope':'Directional variational diagnosis within saved RI-JK scalar-X2C UHF references; not a physical ground-state proof or a production-state replacement.','states':[]}
for z,level in targets:
    start=time.time();path=next(source.glob(f'{z:03d}-*-v{level}z.json'))
    mf,raw=a.load_reference(path)
    tensor=np.array(raw['response']['relaxed_tensor'])
    eigenvalues,directions=np.linalg.eigh(tensor)
    direction=directions[:,0]
    dip=np.einsum('x,xij->ij',direction,a.dipole_matrix(mf))
    ov=mf.get_ovlp();ls=a.angular_labels(mf.mol)
    parity_ao=ov*((-1.)**ls)[None,:]
    cv=[];ci=[];gaps=[];couplings=[];masks=[];parity_min=[]
    for c,e,o in zip(mf.mo_coeff,mf.mo_energy,mf.mo_occ):
        occ=o>0;cv.append(c[:,~occ]);ci.append(c[:,occ])
        gaps.append(e[~occ,None]-e[occ])
        parity=np.einsum('pi,pq,qi->i',c,parity_ao,c,optimize=True)
        parity_min.append(float(np.min(np.abs(parity))))
        assert parity_min[-1]>.99999
        masks.append(parity[~occ,None]*parity[occ]<0)
        couplings.append((cv[-1].T@dip@ci[-1])*masks[-1])
    split=gaps[0].size;shapes=[g.shape for g in gaps]
    def unpack(vector):return [vector[:split].reshape(shapes[0]),vector[split:].reshape(shapes[1])]
    gap=np.concatenate([g.ravel() for g in gaps]);b=np.concatenate([v.ravel() for v in couplings])
    induced=mf.gen_response(hermi=1)
    def action(vector):
        rotation=unpack(vector);densities=[]
        for v,i,u,mask in zip(cv,ci,rotation,masks):
            q=v@(u*mask)@i.T;densities.append(q+q.T)
        potential=induced(np.array(densities))
        return gap*vector+np.concatenate([((v.T@w@i)*mask).ravel() for v,w,i,mask in zip(cv,potential,ci,masks)])
    operator=LinearOperator((len(gap),len(gap)),matvec=action)
    preconditioner=LinearOperator(operator.shape,matvec=lambda v:v/(np.abs(gap)+.2))
    print('START',path.name,'variables',len(gap),flush=True)
    cache=out.with_name(path.stem+'-response-direction.npz')
    cached=cache.exists()
    if cached:
        with np.load(cache,allow_pickle=False) as saved:
            assert str(saved['source_sha256'])==hashlib.sha256(path.read_bytes()).hexdigest()
            u=saved['u'];info=int(saved['minres_info'])
    else:
        u,info=minres(operator,-b,M=preconditioner,rtol=1e-12,maxiter=200)
    residual=float(np.max(np.abs(action(u)+b)/(np.abs(gap)+.2)))
    refine=None
    if residual>2e-9:
        u,refine=gmres(operator,-b,x0=u,M=preconditioner,rtol=1e-11,atol=1e-12,restart=35,maxiter=8)
        residual=float(np.max(np.abs(action(u)+b)/(np.abs(gap)+.2)))
    np.savez_compressed(cache,u=u,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),minres_info=info)
    q=u/np.linalg.norm(u);hq=action(q)
    curvature=2*float(q@hq)
    chi=-2*float(b@u)
    assert residual<2e-8,(path.name,residual)
    assert abs(chi-eigenvalues[0])<max(2e-4,abs(chi)*2e-6),(chi,eigenvalues[0])
    assert curvature<0
    h0=mf.get_hcore();e0=float(mf.energy_tot(dm=mf.make_rdm1(),h1e=h0))
    trials=[]
    if args.stable_energy_differences:
        # UHF energy is quadratic in the density: dE=Tr(F0 dD)+Tr(dD V[dD])/2.
        # An SVD of each occupied/virtual rotation gives exact sin/cos density
        # changes, avoiding subtraction of nearly equal occupied projectors.
        f0=h0+mf.get_veff(mf.mol,mf.make_rdm1())
        stable_modes=[]
        for v,i,rotation,f in zip(cv,ci,unpack(q),f0):
            left,singular,right_t=np.linalg.svd(rotation,full_matrices=False)
            vr=v@left;ir=i@right_t.T
            diagonal_gap=np.einsum('pi,pq,qi->i',vr,f,vr,optimize=True)-np.einsum('pi,pq,qi->i',ir,f,ir,optimize=True)
            cross=np.einsum('pi,pq,qi->i',vr,f,ir,optimize=True)
            stable_modes.append((vr,ir,singular,diagonal_gap,cross))
        def stable_changes(step):
            evens=[];odds=[];linear_even=0.;linear_odd=0.
            for vr,ir,singular,diagonal_gap,cross in stable_modes:
                sine=np.sin(step*singular);cosine=np.cos(step*singular)
                evens.append((vr*sine**2)@vr.T-(ir*sine**2)@ir.T)
                off=(vr*(sine*cosine))@ir.T
                odds.append(off+off.T)
                linear_even+=float(diagonal_gap@(sine**2))
                linear_odd+=2*float(cross@(sine*cosine))
            even=np.array(evens);odd=np.array(odds)
            changes=[]
            for sign in [1,-1]:
                delta=even+sign*odd
                change=linear_even+sign*linear_odd+.5*float(np.einsum('sij,sji->',delta,induced(delta)))
                changes.append(change)
            return changes
    # Compare the derivative at successively halved rotations; retain nonlinear coarse-step differences.
    for step in [.01,.005,.0025]:
        energies=[];counts=[]
        for sign in [1,-1]:
            densities=[]
            for c,o,rotation in zip(mf.mo_coeff,mf.mo_occ,unpack(q)):
                oi=np.flatnonzero(o>0);vi=np.flatnonzero(o==0)
                k=np.zeros((len(o),len(o)))
                k[np.ix_(vi,oi)]=sign*step*rotation
                k[np.ix_(oi,vi)]=-sign*step*rotation.T
                rotated=c@expm(k)
                densities.append((rotated*o)@rotated.T)
            densities=np.array(densities)
            energies.append(float(mf.energy_tot(dm=densities,h1e=h0)))
            counts.append(float(np.einsum('sij,ji->',densities,ov)))
        finite=(sum(energies)-2*e0)/step**2
        trials.append({'rotation_norm':step,'energies':energies,'symmetric_energy_change':sum(energies)/2-e0,'finite_curvature':finite,'relative_curvature_difference':abs(finite/curvature-1),'counts':counts})
        if args.stable_energy_differences:
            changes=stable_changes(step)
            stable_finite=sum(changes)/step**2
            trials[-1].update({'stable_energy_changes':changes,
                              'stable_finite_curvature':stable_finite,
                              'stable_relative_curvature_difference':abs(stable_finite/curvature-1),
                              'maximum_total_subtraction_difference':max(abs(change-(energy-e0)) for change,energy in zip(changes,energies))})
            assert stable_finite<0 and max(changes)<0
        assert finite<0
        assert max(abs(n-z) for n in counts)<2e-7
    comparison='stable_relative_curvature_difference' if args.stable_energy_differences else 'relative_curvature_difference'
    # Save unsuccessful comparisons as evidence too; do not certify them.
    evidence={'state':path.name,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'curvature':curvature,'chi':chi,'residual':residual,'finite_rotations':trials,
              'comparison':comparison,'passed':trials[-1][comparison]<.02}
    out.with_name(path.stem+'-rotation-comparison.json').write_text(json.dumps(evidence,indent=2)+'\n')
    assert trials[-1][comparison]<.02,(path.name,comparison,trials[-1][comparison])
    entry={'response_vector_cached':cached,'finest_rotation_curvature_agrees_within_two_percent':True,'z':z,'symbol':raw['symbol'],'level':level,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_state_sha256':hashlib.sha256(path.with_name(path.stem+'-state.npz').read_bytes()).hexdigest(),'response_eigenvalues':eigenvalues.tolist(),'tested_direction':direction.tolist(),'minimum_orbital_parity_magnitude':min(parity_min),'recomputed_directional_response':chi,'scaled_response_residual':residual,'minres_info':int(info),'gmres_info':None if refine is None else int(refine),'normalized_direction_energy_curvature':curvature,'baseline_energy_difference':e0-raw['energy'],'finite_rotations':trials,'conclusion':'Negative finite energy curvature: this saved stationary determinant is a saddle along a dipole-coupled orbital rotation in its RI-JK functional. The negative response is not a physical ground-state polarizability.','seconds':time.time()-start}
    entry['finite_curvature_comparison']=comparison
    if args.stable_energy_differences:
        entry['energy_difference_method']='Exact RI-JK UHF quadratic density increment Tr(F0 dD)+Tr(dD V[dD])/2; SVD sin/cos occupied/virtual rotations. Original total-energy subtractions are retained separately.'
    result['states'].append(entry)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(entry),flush=True)
