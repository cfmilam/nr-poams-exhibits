#!/usr/bin/env python3
"""Reproduce the helium atomic preparation and hydrogen response controls.

Run: python3 scripts/whole_atom_build.py (requires NumPy).
No network, fitted atomic data, named Gaussian basis or external solver.
Reference values are validation data, never arguments to the energy solve.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import whole_atom_helium as solver

HERE=Path(__file__).resolve().parent
CHECKS=[]


def state_normalization():
    return dict(
        applies_to='helium_infinite and helium4_finite; hydrogen exports spectral response only',
        coordinates='Physical lengths divided by a_*. rho_i=zeta*r_i; s,t,u in polynomial terms are rho1+rho2, rho1-rho2, and |rho1-rho2|.',
        ground='psi = zeta^3/(4*pi*a_*^3) * exp[-(rho1+rho2)] * sum c_ljn s^l t^j u^n. The full six-dimensional norm is one.',
        coefficient_metric='The solver suppresses the common 16*pi^2 in the six-dimensional angular/perimetric measure. Exported scalar c satisfies c^T S c=1 in that suppressed measure.',
        response='delta_psi/F = -zeta^3/(4*pi*a_*^3) * exp[-(rho1+rho2)] * sum q_ljn s^l t^j u^n (rho1_vector + sign*rho2_vector)_z. F is dimensionless; q are response.coefficients.',
        intrinsic='Only the normalized spatial amplitude is exported; the intrinsic singlet factor is specified separately.',
        spectral='gaps are in E_*; strengths are squared one-Cartesian matrix elements of D/a_*. They are invariant to basis scaling and directly usable in the response/C6 formulas.')


def check(name, condition, actual=None, tolerance=None):
    record=dict(name=name,passed=bool(condition))
    if actual is not None: record['actual']=float(actual)
    if tolerance is not None: record['tolerance']=tolerance
    CHECKS.append(record)
    if not condition: raise AssertionError(record)


def compare(name,a,b,tol):
    check(name,abs(a-b)<tol,abs(a-b),tol)


def brief(result):
    row={k:result[k] for k in ('degree','quadrature','energy','zeta','support','removal_interval')}
    row.update(polarizability=result['response']['polarizability'],
               C6=result['response']['C6_identical'],response_support=result['response']['support'])
    return row


def frequency_C6(a,b):
    t,w=np.polynomial.legendre.leggauss(120)
    t=(t+1)/2; w=w/2; om=t/(1-t)
    def response(r):
        e=np.asarray(r['gaps']); s=np.asarray(r['strengths'])
        return 2*np.sum(e[:,None]*s[:,None]/(e[:,None]**2+om[None,:]**2),axis=0)
    return float(3/math.pi*np.sum(w*response(a)*response(b)/(1-t)**2))


def main():
    shared=json.loads((HERE/'coulomb_calibrated_inputs.json').read_text())['inputs']
    val=lambda k:float(shared[k]['value'])
    h,c,alpha,m,e=map(val,('h','c','alpha','m_light','e_SI'))
    M=6.6446573450e-27
    eta=m/M
    a=h/(2*math.pi*m*c*alpha); E=m*c*c*alpha*alpha; rate=E/h
    units=dict(length_m=a,energy_J=E,energy_eV=E/e,rate_Hz=rate,
               force_N=E/a,compliance_m2_per_J=a*a/E,C6_J_m6=E*a**6)
    baseline=solver.solve(0,18,0,response_degree=2)
    baseline_finite=solver.solve(0,18,eta,response_degree=2)
    compare('analytic one-shape energy',baseline['energy'],-(27/16)**2,1e-10)
    compare('analytic one-shape scale',baseline['zeta'],27/16,1e-6)
    compare('analytic uncorrelated mean distance',baseline['moments']['mean_r'],8/9,1e-6)
    compare('analytic uncorrelated separation',baseline['moments']['mean_r12'],35/27,1e-6)

    infinity=[]; finite=[]
    for p in (2,4,6,8,10,12):
        q=max(20,p+14)
        infinity.append(solver.solve(p,q,0,response_degree=p))
        finite.append(solver.solve(p,q,eta,response_degree=p))
        print('degree',p,'finite energy',finite[-1]['energy'],
              'response',finite[-1]['response']['polarizability'],flush=True)
    inf,fin=infinity[-1],finite[-1]
    # Degree, quadrature and retained-support perturbations are independent
    # numerical checks, not new physics observations or parameter fits.
    variants={
        'quadrature32':solver.solve(12,32,eta,response_degree=12),
        'cutoff1e-12':solver.solve(12,26,eta,1e-12,12),
        'cutoff1e-14':solver.solve(12,26,eta,1e-14,12),
        'degree14':solver.solve(14,28,eta,response_degree=14),
    }
    for name,r in variants.items():
        compare(name+': energy stability',r['energy'],fin['energy'],1e-8)
        compare(name+': static response stability',r['response']['polarizability'],fin['response']['polarizability'],3e-6)
        compare(name+': dispersion stability',r['response']['C6_identical'],fin['response']['C6_identical'],3e-6)
        compare(name+': envelope mean radius stability',r['moments']['mean_r'],fin['moments']['mean_r'],2e-6)
    for label,rows in [('infinite',infinity),('finite',finite)]:
        check(label+': decreasing calculated energies',all(rows[i+1]['energy']<rows[i]['energy']+1e-10 for i in range(len(rows)-1)))
        r=rows[-1]; pol=r['response']
        compare(label+': energy sum',sum(r[k] for k in ('light_kinetic','eye_recoil','eye_attraction','mobile_pair_repulsion')),r['energy'],1e-8)
        compare(label+': virial',r['virial'],0,2e-7)
        compare(label+': normalization',r['normalization'],1,1e-9)
        compare(label+': dipole strength completeness',pol['dipole_strength_sum'],r['moments']['dipole_z_squared'],1e-8)
        compare(label+': energy-weighted dipole sum',pol['energy_weighted_sum'],pol['commutator_reference'],2e-7)
        compare(label+': independent frequency-integral C6',frequency_C6(pol,pol),pol['C6_identical'],2e-9)
        check(label+': all response denominators positive',min(pol['gaps'])>0)
        check(label+': positive recoil and correlation variances',r['eye_recoil']>=0 and r['moments']['relative_one_entry_L_squared_over_hbar_squared']>0)
    refs=dict(infinite_energy=-2.9037243770341195983,
              infinite_static_response=1.383192174455,
              finite_nonrelativistic_static_response=1.383809986408,
              corrected_physical_static_response=1.38376078,
              helium_first_removal_eV=24.587389011,
              H_H_C6=6.499026705405839,
              infinite_He_He_C6=1.46097783768,
              infinite_H_He_C6=2.82134391528)
    compare('independent infinite-eye energy',inf['energy'],refs['infinite_energy'],1e-8)
    compare('independent infinite-eye response',inf['response']['polarizability'],refs['infinite_static_response'],3e-6)
    compare('independent finite-eye response',fin['response']['polarizability'],refs['finite_nonrelativistic_static_response'],3e-6)
    compare('independent infinite-eye dispersion',inf['response']['C6_identical'],refs['infinite_He_He_C6'],3e-6)
    hydrogen={str(p):solver.hydrogen_response(p) for p in (0,1,4,8,10,12)}
    H=hydrogen['10']
    compare('hydrogen occupied-only static response',H['occupied_only_response'],0,1e-15)
    compare('hydrogen one directional radial function',hydrogen['0']['polarizability'],4,1e-12)
    compare('hydrogen exact static response with two radial functions',hydrogen['1']['polarizability'],4.5,1e-12)
    compare('hydrogen converged static response',H['polarizability'],4.5,1e-9)
    compare('hydrogen total dipole strength',H['dipole_strength_sum'],1,1e-9)
    compare('hydrogen energy-weighted sum',H['energy_weighted_sum'],.5,1e-9)
    compare('hydrogen dispersion benchmark',H['C6_identical'],refs['H_H_C6'],1e-8)
    compare('mixed infinite-eye H/He dispersion',solver.dispersion(H,inf['response']),refs['infinite_H_He_C6'],4e-6)
    compare('hydrogen dispersion degree convergence',H['C6_identical'],hydrogen['12']['C6_identical'],1e-8)
    # A two-sector finite-field calculation independently reads the curvature
    # of the combined matrix, instead of repeating the inverse-matrix formula.
    gaps=np.asarray(fin['response']['gaps']); d=np.sqrt(fin['response']['strengths'])
    block=np.diag(np.r_[0.,gaps]); F=1e-4
    block[0,1:]=F*d;block[1:,0]=F*d
    curvature=-2*np.linalg.eigvalsh(block)[0]/F**2
    compare('finite-field energy curvature',curvature,fin['response']['polarizability'],1e-6)

    radii=np.linspace(0,8,401)
    probability=solver.radial_profile(fin,radii)
    compare('independent radial marginal normalization',np.trapezoid(probability,radii),1,5e-5)
    compare('independent radial marginal mean radius',np.trapezoid(probability*radii,radii),fin['moments']['mean_r'],5e-5)
    compare('independent radial marginal second moment',np.trapezoid(probability*radii*radii,radii),fin['moments']['mean_r2'],1e-4)
    check('radial distribution nonnegative',min(probability)>=0)
    cumulative=np.r_[0,np.cumsum((probability[1:]+probability[:-1])*np.diff(radii)/2)]
    radius90=float(np.interp(.9,cumulative,radii))
    energy_rows=[]
    for key,title in [('light_kinetic','Mobile kinetic sum'),('eye_recoil','Compact-component translational recoil'),
                      ('eye_attraction','Joint compact/mobile attraction'),('mobile_pair_repulsion','Joint mobile/mobile repulsion'),
                      ('energy','Total external binding account'),('removal_interval','First-removal interval')]:
        value=fin[key]
        energy_rows.append(dict(key=key,label=title,value=value,rate_PHz=value*rate/1e15,energy_J=value*E,energy_eV=value*E/e))
    ion=fin['removal_interval']*E/e
    results=dict(
        schema='poams.whole_atom.helium.v1',date='2026-09-15',
        model='Effective nonrelativistic inverse-distance helium operator; measured alpha; finite point-eye recoil',
        status='Correlated envelope and response solved within the declared model. Intrinsic compact-core structure not derived.',
        input_policy='No helium energy, radius, polarizability or dispersion coefficient is fitted. Species count and symmetry sector are specified, not derived from isotope counts.',
        inputs=dict(alpha=alpha,h=h,c=c,m_light_kg=m,M_eye_kg=M,M_eye_uncertainty_kg=2.1e-36,eta=eta,Z=2,mobile_count=2,
                    compact_I=0,compact_I_status='measured state assignment',intrinsic_pair='singlet sector specified',
                    source='https://physics.nist.gov/cuu/Constants/Table/allascii.txt'),
        units=units,references=refs,
        sources=[
            dict(title='NIST ASD helium ionization value',url='https://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=2',role='observable comparison, not input'),
            dict(title='Schwartz 2006 helium energy',url='https://arxiv.org/abs/math-ph/0605018',role='same-operator numerical benchmark'),
            dict(title='Puchalski, Jentschura and Mohr 2011',url='https://www.nist.gov/document/blackbody-radiationpdf',role='infinite-eye response benchmark'),
            dict(title='Puchalski et al. 2020 helium response',url='https://arxiv.org/abs/1912.12242',role='finite-mass and corrected-response comparisons'),
            dict(title='Yan, Babb, Dalgarno and Drake 1996',url='https://arxiv.org/abs/atom-ph/9607002',role='H/H, H/He and He/He dispersion benchmarks'),
            dict(title='NIST helium isotope and ground state',url='https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable1.htm',role='compact I=0 state input and angular-state correspondence')],
        convergence=dict(infinite=[brief(r) for r in infinity],finite=[brief(r) for r in finite],variants={k:brief(v) for k,v in variants.items()}),
        uncorrelated_baseline=dict(energy=baseline_finite['energy'],mean_radius_pm=baseline_finite['moments']['mean_r']*a*1e12,
                                   first_removal_eV=baseline_finite['removal_interval']*E/e,zeta=baseline_finite['zeta']),
        energy_account=energy_rows,
        envelope=dict(**fin['moments'],mean_pm=fin['moments']['mean_r']*a*1e12,rms_pm=fin['moments']['rms_r']*a*1e12,
                      mean_separation_pm=fin['moments']['mean_r12']*a*1e12,radius90_pm=radius90*a*1e12,
                      radius90_status='numerical radial quantile, not a hard atomic boundary',
                      radial_profile=[dict(radius_pm=float(r*a*1e12),probability_per_pm=float(p/(a*1e12))) for r,p in zip(radii,probability)]),
        response=dict(finite=fin['response']['polarizability'],infinite=inf['response']['polarizability'],
                      occupied_only=0.,numerical_display_tolerance=3e-6,
                      physical_reference_difference_ppm=(fin['response']['polarizability']/refs['corrected_physical_static_response']-1)*1e6,
                      compliance_m2_per_J=fin['response']['polarizability']*a*a/E,
                      commutator_sum=fin['response']['energy_weighted_sum'],commutator_reference=1+2*eta),
        dispersion=dict(H_H=H['C6_identical'],infinite_H_He=solver.dispersion(H,inf['response']),
                        infinite_He_He=inf['response']['C6_identical'],finite_He4_He4=fin['response']['C6_identical'],
                        units='E_* a_*^6',scope='Non-overlapping, nonretarded leading dipole response only; no bond curve, short-range repulsion or dimer binding predicted.'),
        comparisons=dict(first_removal_eV=ion,first_removal_reference_eV=refs['helium_first_removal_eV'],
                         first_removal_residual_eV=ion-refs['helium_first_removal_eV'],
                         first_removal_residual_ppm=(ion/refs['helium_first_removal_eV']-1)*1e6,
                         infinite_energy_error=inf['energy']-refs['infinite_energy']),
        hydrogen_controls={k:{n:v for n,v in r.items() if n not in ('gaps','strengths')} for k,r in hydrogen.items()},
        checks=CHECKS,
        numerical_precision='Display about 5-6 significant digits for envelope response. Degree/support checks indicate 1e-8 E_* for energy and 3e-6 absolute for response/C6; these are convergence indicators, not rigorous error bars or physical uncertainty.',
        not_derived=['compact-core shape, radius, thickness and intrinsic binding','regional core/envelope energy or count-to-mass law',
                     'carbon/nitrogen/oxygen many-entry response','molecular short-range exchange and bond geometry','propulsion or whole-device support response'],
        molecular_use=['exact hydrogen directional-response gate','two-entry joint-tensor and recoil benchmark',
                       'correlated helium envelope moments','complete dipole response and energy-weighted sums',
                       'leading long-range H/H, H/He and He/He C6 controls'],
        next_gate='Apply corresponding complete response spaces to the existing C/N/O atom generators; retain water representation/gradient gates before any long geometry search.',
    )
    handoff=dict(schema='poams.atomic_response_handoff.v1',units=units,inputs=results['inputs'],
                 state_normalization=state_normalization(),
                 use='Benchmarks and long-range response; not a plug-compatible replacement for C/N/O four-index tensors.',
                 hydrogen=H,helium_infinite=inf,helium4_finite=fin,
                 warning='Do not interpret finite response pseudostates as a converged measured excitation spectrum. Do not add C6 to a correlated molecular operator already containing the same response.',
                 checks_file='whole_atom_results.json')
    (HERE/'whole_atom_results.json').write_text(json.dumps(results,indent=2)+'\n')
    display={k:results[k] for k in ('date','status','units','energy_account','envelope','response','dispersion','comparisons','uncorrelated_baseline')}
    (HERE/'whole_atom_data.js').write_text('// Generated by whole_atom_build.py; do not hand-edit.\nwindow.POAMS_WHOLE_ATOM = '+json.dumps(display,separators=(',',':'))+';\n')
    (HERE/'whole_atom_handoff.json').write_text(json.dumps(handoff,indent=2)+'\n')
    with (HERE/'whole_atom_helium_radial.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['radius_pm','one_entry_probability_per_pm','cumulative_probability'])
        writer.writerows(zip(radii*a*1e12,probability/(a*1e12),cumulative))
    print(json.dumps(dict(checks=len(CHECKS),envelope_pm=results['envelope']['mean_pm'],response=results['response'],
                          dispersion=results['dispersion'],comparisons=results['comparisons']),indent=2))


if __name__=='__main__': main()
