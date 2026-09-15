#!/usr/bin/env python3
"""Source-grounded circulation composition and its selection boundary.

Osborne & Pope (2007), chapter 7, postulate an orbital kinetic-energy norm for
ONE axially symmetric spinning body. This is not the exhibit's later bilinear
-2G J1.J2/(c^2 r^3) ansatz. Derive the source ratio, then give two explicit
multi-circulation continuations that agree on the source's single-body case.
Neither continuation is selected as the superconducting device's physics.

All positions are relational chart labels. The tensor-kernel examples below
test symmetry, reciprocal reactions and composition, not a field or mediator.
No experimental force is fitted. NumPy is the only third-party dependency.
"""
from decimal import Decimal, localcontext
import json
import math

import numpy as np

CHECKS = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def near(a, b, rtol=1e-10, atol=1e-12):
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))


def source_ratio(epsilon, cosine):
    """Gcal/G = vs/v0 = r0/rs, with Ko>0, Ks>=0 and no zero-energy cusp."""
    return (1 + 2*epsilon*cosine + epsilon*epsilon)**.25


def single_excess(epsilon, cosine):
    # Avoid subtracting two near-unit floating-point values.
    return math.expm1(.25*math.log1p(2*epsilon*cosine + epsilon*epsilon))


def paired_excess(epsilon, cosine):
    """Delta N/(total rotor mass * GM_E/R^2), entrywise continuation only.

    Decimal protects cancellation of two O(epsilon) terms to O(epsilon^2).
    These illustrative inputs are dimensionless, not apparatus measurements.
    """
    with localcontext() as ctx:
        ctx.prec = 65
        e, u = Decimal(str(epsilon)), Decimal(str(cosine))
        plus = (1+e*e+2*e*u).sqrt().sqrt()
        minus = (1+e*e-2*e*u).sqrt().sqrt()
        return float((plus+minus)/2-1)


def local_load(masses, ko, ks, axes, orbital_axis, g0=1.):
    """Sum separately evaluated source-body contributions; added premise."""
    return g0*sum(m*single_excess(s/o, float(n@orbital_axis))
                  for m, o, s, n in zip(masses, ko, ks, axes))


def collective_load(masses, ko, ks, axes, orbital_axis, g0=1.):
    """Combine directional energy entries before taking norm; alternative.

    Not an identification of counter-rotation with a spinless rigid body:
    all positive stored kinetic energies remain in the retained state list.
    """
    orbital = np.sum(ko)
    w = np.sum(ks[:, None]*axes, axis=0)/orbital
    return g0*np.sum(masses)*math.expm1(.25*math.log1p(2*float(w@orbital_axis)+float(w@w)))


def j_total(positions, momenta, spins):
    return np.sum(np.cross(positions, momenta)+spins, axis=0)


def internal_book(masses, positions, momenta, spins):
    mass = np.sum(masses)
    center = np.sum(masses[:, None]*positions, axis=0)/mass
    momentum = np.sum(momenta, axis=0)
    local_p = momenta-masses[:, None]*momentum/mass
    inner = np.sum(np.cross(positions-center, local_p)+spins, axis=0)
    return center, momentum, inner


def tensor_energy(rvec, j, k, aa=-2., bb=0., coefficient=1.):
    """Restricted reciprocal, parity-even, bilinear, no-new-length family."""
    r = np.sqrt(rvec@rvec)
    n = rvec/r
    return coefficient*(aa*(j@k)+bb*(j@n)*(k@n))/r**3


def tensor_force(rvec, j, k, aa=-2., bb=0., coefficient=1.):
    r = np.sqrt(rvec@rvec)
    n = rvec/r
    jj, kk = j@n, k@n
    return coefficient*(3*aa*(j@k)*n+5*bb*jj*kk*n-bb*(kk*j+jj*k))/r**4


def torques(rvec, j, k, aa, bb, coefficient=1.):
    r = np.sqrt(rvec@rvec)
    n = rvec/r
    return (-coefficient*np.cross(j, aa*k+bb*(k@n)*n)/r**3,
            -coefficient*np.cross(k, aa*j+bb*(j@n)*n)/r**3)


def cross_cost(x, js, y, ks, aa, bb):
    return sum(tensor_energy(xi-yi, ji, ki, aa, bb)
               for xi, ji in zip(x, js) for yi, ki in zip(y, ks))


def derivative(fn):
    return np.imag(fn(1e-25j))/1e-25


def main():
    # Reconstruct source equations 7.11-7.21 independently of the ratio function.
    mass, primary_mass, grav, action = 2., 100., .7, 10.
    v0 = grav*mass*primary_mass/action
    r0 = action/(mass*v0)
    ko = mass*v0*v0/2
    for e, cosine in [(0., .3), (.2, 1.), (.2, -1.), (.2, 0.), (.37, .6)]:
        ks = e*ko
        total_k = math.sqrt(ko*ko+ks*ks+2*ko*ks*cosine)
        vs = math.sqrt(2*total_k/mass)
        rs = action/(mass*vs)
        gc = rs*vs*vs/primary_mass
        h = source_ratio(e, cosine)
        check(f'source circular equations reconstruct fourth-root law {e, cosine}',
              near([vs/v0, r0/rs, gc/grav], [h]*3))
    check('parallel source is sqrt(1+epsilon)', near(source_ratio(.2, 1), math.sqrt(1.2)))
    check('opposed source is sqrt(1-epsilon) below the cusp', near(source_ratio(.2, -1), math.sqrt(.8)))
    check('perpendicular source has a nonzero even term', paired_excess(.2, 0)>0)
    check('single-source small-spin term is cosine*epsilon/2',
          near(single_excess(1e-7, .4)/1e-7, .2, rtol=1e-6))
    for cosine in [0., .3, .7, 1.]:
        e=.001
        leading = (2-3*cosine*cosine)*e*e/8
        residual = paired_excess(e, cosine)-leading
        half_residual = paired_excess(e/2, cosine)-(2-3*cosine*cosine)*(e/2)**2/8
        check(f'paired source second-order term with fourth-order remainder {cosine}',
              near(residual/half_residual,16.,rtol=2e-5))
    check('matched source pair is even in orientation', near(paired_excess(.2,.7), paired_excess(.2,-.7)))
    check('parallel source pair has negative second-order excess', paired_excess(.2,1)<0)
    check('equal-energy reversal gives no source-model switch signal',
          paired_excess(.2,.7)-paired_excess(.2,-.7)==0)
    check('changed stored kinetic energy can change paired source load',
          paired_excess(.2,0)-paired_excess(.1,0)>0)
    check('symmetric AM-sign history need not average even load to zero',
          (paired_excess(.2,.7)+paired_excess(.2,-.7))/2!=0)
    # Newtonian-equivalent support term: the held-rate contribution cancels
    # in spin-state subtraction, not in the absolute baseline weight.
    R, held_rate = 4., .03
    gp = R**3*held_rate**2/primary_mass
    gc = grav*source_ratio(.2,.3)
    n0 = mass*primary_mass*(grav-gp)/R**2
    ns = mass*primary_mass*(gc-gp)/R**2
    check('source support subtraction keeps the same held constraint',
          near(ns-n0, mass*grav*primary_mass/R**2*single_excess(.2,.3)))

    masses=np.array([1.,1.]); orbital=np.array([2.,2.]); stored=np.array([.4,.4])
    z=np.array([0.,0.,1.]); axes=np.array([[1.,0.,0.],[-1.,0.,0.]])
    nl=local_load(masses,orbital,stored,axes,z)
    nc=collective_load(masses,orbital,stored,axes,z)
    check('two compositions agree for every named single rotor',
          all(near(local_load(masses[:1],orbital[:1],stored[:1],axes[:1],n),
                   collective_load(masses[:1],orbital[:1],stored[:1],axes[:1],n))
              for n in [z,np.array([1.,0.,0.]),np.array([.6,0.,.8])]))
    check('source-consistent compositions disagree for opposed rotors', nl>0 and nc==0)
    check('entrywise opposite-rotor sum equals analytic pair formula',near(nl/2,paired_excess(.2,0)))
    same=np.array([[1.,0.,0.],[1.,0.,0.]])
    check('both continuations agree for a homogeneous co-rotating body',
          near(local_load(masses,orbital,stored,same,z),collective_load(masses,orbital,stored,same,z)))
    for fn in [local_load,collective_load]:
        # Split an existing homogeneous collection, retaining its energy/axis
        # metadata; no fractional physical closure is postulated.
        ms=np.array([.25,.75,1.]); os=np.array([.5,1.5,2.]); ss=np.array([.1,.3,.4])
        nn=np.array([axes[0],axes[0],axes[1]])
        check(fn.__name__+' homogeneous accounting subdivision invariant',
              near(fn(ms,os,ss,nn,z), fn(masses,orbital,stored,axes,z)))
        check(fn.__name__+' depends on energy ratios not common energy units',
              near(fn(masses,orbital*5,stored*5,axes,z),fn(masses,orbital,stored,axes,z)))
        check(fn.__name__+' common rotation of both senses unchanged',
              near(fn(masses,orbital,stored,-axes,-z),fn(masses,orbital,stored,axes,z)))

    # Actual AM is additive independently of any proposed energy law.
    masses=np.array([2.,3.,5.]); pos=np.array([[1.,2.,0.],[-2.,1.,1.],[1.,0.,3.]])
    mom=np.array([[2.,0.,1.],[-1.,3.,0.],[0.,2.,-1.]])
    spin=np.array([[0.,0.,2.],[1.,0.,0.],[0.,-1.,0.]])
    ctr,p,j=internal_book(masses,pos,mom,spin)
    check('COM orbital plus internal AM equals complete finite sum',near(j_total(pos,mom,spin),np.cross(ctr,p)+j))
    groups=[np.array([0,1]),np.array([2])]
    grouped=[]
    for ix in groups:
        cr,pr,jr=internal_book(masses[ix],pos[ix],mom[ix],spin[ix])
        grouped.append(np.cross(cr,pr)+jr)
    check('grouping retains local orbital and intrinsic entries exactly',near(sum(grouped),j_total(pos,mom,spin)))
    check('zero net opposed AM does not imply zero stored kinetic energy',
          np.sum(np.array([[0.,0.,2.],[0.,0.,-2.]]),axis=0).tolist()==[0.,0.,0.] and 2*(2**2/(2*3))>0)
    # J = -dB_rot/dOmega at zero rotation does not fix the kinetic decomposition.
    s, inertias = 2., [3.,7.]
    kinetic=[s*s/(2*i) for i in inertias]
    total=10.
    check('identical J and total energy admit different kinetic partitions',
          kinetic[0]!=kinetic[1] and all(near(k+(total-k),total) for k in kinetic))

    r=np.array([1.4,-.3,.7]); j=np.array([.2,.7,-.5]); k=np.array([-.6,.4,.8])
    rot,_=np.linalg.qr(np.array([[1.,2.,3.],[4.,0.,1.],[2.,1.,5.]]))
    for aa,bb in [(-2.,0.),(-2.,1.)]:
        f=tensor_force(r,j,k,aa,bb)
        numeric=-np.array([derivative(lambda t: tensor_energy(r+t*e,j,k,aa,bb)) for e in np.eye(3)])
        check(f'tensor family {aa,bb}: analytic reaction agrees with energy derivative',near(f,numeric))
        check(f'tensor family {aa,bb}: endpoint exchange gives equal/opposite reaction',near(tensor_force(-r,k,j,aa,bb),-f))
        t1,t2=torques(r,j,k,aa,bb)
        check(f'tensor family {aa,bb}: orbital and local torques sum to zero',near(np.cross(r,f)+t1+t2,0.))
        check(f'tensor family {aa,bb}: common rotation invariant',near(tensor_energy(rot@r,rot@j,rot@k,aa,bb),tensor_energy(r,j,k,aa,bb)))
        check(f'tensor family {aa,bb}: parity and time reversal invariant',
              near(tensor_energy(-r,j,k,aa,bb),tensor_energy(r,-j,-k,aa,bb)))
        check(f'tensor family {aa,bb}: distance scale is r^-3 under stated restrictions',near(tensor_energy(2*r,j,k,aa,bb),tensor_energy(r,j,k,aa,bb)/8))
    check('two admissible tensor kernels are different physical laws',
          not near(tensor_energy(r,j,k,-2,0),tensor_energy(r,j,k,-2,1)))
    n=np.array([0.,0.,1.]); transverse=np.array([1.,0.,0.])
    check('one transverse calibration cannot choose longitudinal coefficient',
          tensor_energy(n,transverse,transverse,-2,0)==tensor_energy(n,transverse,transverse,-2,1))
    # r^-3 follows dimensionally only if G/c^2 and two AM entries are the only
    # dimensional building blocks. Dimensions below are (mass,length,time).
    check('G/c^2 * J^2 / r^3 has energy dimensions',
          np.array_equal(np.array([-1,1,0])+2*np.array([1,2,-1])-np.array([0,3,0]),[1,2,-2]))
    x=np.array([[0.,0.,0.],[.1,.2,0.]]); js=np.array([[1.,0.,1.],[0.,1.,-.2]])
    y=np.array([[0.,0.,3.],[.3,0.,4.]]); ks=np.array([[.2,.1,1.],[.4,-.5,.2]])
    for aa,bb in [(-2.,0.),(-2.,1.)]:
        val=cross_cost(x,js,y,ks,aa,bb)
        grouped=sum(cross_cost(x[ix:ix+1],js[ix:ix+1],y,ks,aa,bb) for ix in range(2))
        check(f'finite tensor sum {aa,bb}: exact grouping without missing cross terms',near(val,grouped))
        xs=np.vstack([x[0],x[0],x[1]]); jss=np.vstack([js[0]*.4,js[0]*.6,js[1]])
        check(f'finite tensor sum {aa,bb}: co-located record subdivision invariant',near(val,cross_cost(xs,jss,y,ks,aa,bb)))
    values=np.array([[0.,2.,3.],[2.,0.,5.],[3.,5.,0.]])
    check('co-signed endpoints do not by counting alone fix the physical coefficient',
          np.sum(values)/2==np.sum(np.triu(values,1))==10.)

    examples={str(angle):{'cosine':math.cos(math.radians(angle)),
              'entrywise_delta_N_over_total_mass_g':paired_excess(.2,math.cos(math.radians(angle))),
              'collective_delta_N_over_total_mass_g':0.}
              for angle in [0,35.264389682754654,45,90]}
    result={'status':'conditional source law; multi-circulation and material maps unselected',
      'source':'Osborne and Pope 2007, chapter 7, printed pp.133-135 and eq.7.21',
      'checks_passed':len(CHECKS),'checks':CHECKS,
      'illustrative_dimensionless_epsilon':.2,'examples':examples,
      'perpendicular_epsilon_1e-7':paired_excess(1e-7,0),
      'no_apparatus_force_prediction':True}
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
