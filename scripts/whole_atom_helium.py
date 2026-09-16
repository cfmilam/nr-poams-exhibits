#!/usr/bin/env python3
"""Explicitly correlated two-entry effective atom. Python + NumPy only.

All integrals are computed here; no helium energy or response is a fit input.
Lengths and energies are in calibrated a_* and E_* units. The numerical
operator is the nonrelativistic point-eye helium Hamiltonian, not a derivation
of compact core structure or of the operator from POAMS ontology.
"""
import argparse
import json
import math
import time
from pathlib import Path
import numpy as np


def powers(degree, odd=False):
    # Monomials in (s,t,u); t parity gives exchange parity.
    return [(l, j, n) for total in range(degree + 1)
            for j in range(int(odd), total + 1, 2)
            for l in range(total - j + 1)
            for n in [total - j - l]]


def grid(q):
    a, w = np.polynomial.laguerre.laggauss(q)
    x, y, z = np.meshgrid(a / 2, a / 2, a / 4, indexing='ij')
    wx, wy, wz = np.meshgrid(w, w, w, indexing='ij')
    x, y, z = x.ravel(), y.ravel(), z.ravel()
    r1, r2, u = x+z, y+z, x+y
    weight = (wx*wy*wz).ravel()*r1*r2*u/16
    # Common factor 16*pi^2 cancels all normalized results.
    cos = np.clip((r1*r1+r2*r2-u*u)/(2*r1*r2), -1, 1)
    sin = np.sqrt(np.maximum(0, 1-cos*cos))
    v1 = np.column_stack((r1, np.zeros_like(r1), np.zeros_like(r1)))
    v2 = np.column_stack((r2*cos, r2*sin, np.zeros_like(r2)))
    return dict(r1=r1, r2=r2, u=u, s=r1+r2, t=r1-r2, w=weight,
                cos=cos, v1=v1, v2=v2, h1=v1/r1[:, None],
                h2=v2/r2[:, None], hu=(v1-v2)/u[:, None])


def monomials(g, terms):
    coords = [g[k] for k in ('s', 't', 'u')]
    p, ds, dt, du = [], [], [], []
    def mono(exponents):
        value = np.ones(len(g['s']))
        for x, n in zip(coords, exponents):
            value *= x**n
        return value
    for term in terms:
        value = mono(term)
        p.append(value)
        for k, out in enumerate((ds, dt, du)):
            if term[k]:
                e = list(term)
                e[k] -= 1
                out.append(term[k]*mono(e))
            else:
                out.append(np.zeros_like(value))
    P, S, T, U = [np.array(v).T for v in (p, ds, dt, du)]
    # exp(-s) is factored into Gauss-Laguerre weights.
    return P, S+T-P, S-T-P, U


def gram(a, b, w):
    return a.T @ (w[:, None]*b)


def scalar_matrices(g, degree, Z=2):
    terms = powers(degree)
    P, D1, D2, Du = monomials(g, terms)
    w = g['w']
    S = gram(P, P, w)
    T = np.zeros_like(S)
    R = np.zeros_like(S)
    for k in range(3):
        d1 = D1*g['h1'][:, k, None] + Du*g['hu'][:, k, None]
        d2 = D2*g['h2'][:, k, None] - Du*g['hu'][:, k, None]
        T += (gram(d1, d1, w)+gram(d2, d2, w))/2
        R += gram(d1+d2, d1+d2, w)/2
    Veye = gram(P, P, w*(-Z/g['r1']-Z/g['r2']))
    Vpair = gram(P, P, w/g['u'])
    return dict(terms=terms, P=P, S=S, T=T, recoil=R,
                Veye=Veye, Vpair=Vpair)


def support(S, cutoff):
    diagonal = np.sqrt(np.diag(S))
    if np.any(diagonal <= 0):
        raise ValueError('nonpositive Gram diagonal')
    G = S/diagonal[:, None]/diagonal[None, :]
    e, v = np.linalg.eigh((G+G.T)/2)
    keep = e > cutoff*e[-1]
    X = v[:, keep]/np.sqrt(e[keep])[None, :]/diagonal[:, None]
    return X, dict(dimension=len(e), retained=int(sum(keep)),
                   cutoff=cutoff, smallest_retained_relative=float(e[keep][0]/e[-1]),
                   most_negative_relative=float(min(0, e[0])/e[-1]))


def ground(m, eta=0, cutoff=1e-13, optimize=True, zeta=2.):
    X, info = support(m['S'], cutoff)
    t = X.T @ (m['T']+eta*m['recoil']) @ X
    v = X.T @ (m['Veye']+m['Vpair']) @ X
    def eig(z, vectors=False):
        H = z*z*t+z*v
        H = (H+H.T)/2
        return np.linalg.eigh(H) if vectors else np.linalg.eigvalsh(H)[0]
    if optimize:
        # Energy only. Fixed bracket, no target consulted.
        lo, hi = .8, 3.0
        phi = (math.sqrt(5)-1)/2
        a, b = hi-phi*(hi-lo), lo+phi*(hi-lo)
        ea, eb = eig(a), eig(b)
        for _ in range(55):
            if ea < eb:
                hi, b, eb = b, a, ea
                a = hi-phi*(hi-lo)
                ea = eig(a)
            else:
                lo, a, ea = a, b, eb
                b = lo+phi*(hi-lo)
                eb = eig(b)
        zeta = (lo+hi)/2
    e, c = eig(zeta, True)
    c = X @ c[:, 0]
    if (m['P']@c)[0] < 0:
        c *= -1
    avg = lambda A: float(c @ A @ c)
    Tlight = zeta*zeta*avg(m['T'])
    Teye = zeta*zeta*eta*avg(m['recoil'])
    Ve, Vp = zeta*avg(m['Veye']), zeta*avg(m['Vpair'])
    result = dict(energy=float(e[0]), zeta=zeta, light_kinetic=Tlight,
                  eye_recoil=Teye, eye_attraction=Ve, mobile_pair_repulsion=Vp,
                  virial=2*(Tlight+Teye)+Ve+Vp, support=info,
                  normalization=avg(m['S']), eta=eta,
                  removal_interval=-2/(1+eta)-float(e[0]))
    return result, c


def moments(g, m, c, zeta):
    prob = g['w']*(m['P']@c)**2
    prob /= sum(prob)
    out = {}
    for name, values in [('mean_r',(g['r1']+g['r2'])/2/zeta),
                         ('mean_r2',(g['r1']**2+g['r2']**2)/2/zeta**2),
                         ('mean_r12',g['u']/zeta),
                         ('mean_r12_squared',g['u']**2/zeta**2),
                         ('mean_inverse_r',(1/g['r1']+1/g['r2'])*zeta/2),
                         ('mean_cos_angle',g['cos']),
                         ('dipole_z_squared',np.sum((g['v1']+g['v2'])**2,axis=1)/3/zeta**2)]:
        out[name]=float(prob@values)
    out['rms_r']=math.sqrt(out['mean_r2'])
    _, D1, D2, Du = monomials(g, m['terms'])
    transverse = g['r1']*g['r2']*np.sqrt(1-g['cos']**2)/g['u']
    out['relative_one_entry_L_squared_over_hbar_squared'] = float(
        g['w'] @ ((Du@c)*transverse)**2)
    out['relative_L1_dot_L2_over_hbar_squared'] = -out['relative_one_entry_L_squared_over_hbar_squared']
    out['total_orbital_L_squared_over_hbar_squared'] = 0.0
    return out


def response_matrices(g, degree, Z=2):
    plus, minus = powers(degree), powers(degree, odd=True)
    terms = plus+minus
    signs = np.array([1.]*len(plus)+[-1.]*len(minus))
    P, D1, D2, Du = monomials(g, terms)
    w = g['w']/3
    S = np.zeros((len(terms),len(terms)))
    T, R, V = [np.zeros_like(S) for _ in range(3)]
    potential = -Z/g['r1']-Z/g['r2']+1/g['u']
    rhs_density = np.zeros_like(P)
    for a in range(3):
        vec = g['v1'][:,a,None]+g['v2'][:,a,None]*signs
        f = vec*P
        S += gram(f,f,w)
        V += gram(f,f,w*potential)
        rhs_density += vec*P*(g['v1'][:,a,None]+g['v2'][:,a,None])/3
        for b in range(3):
            d1 = vec*(D1*g['h1'][:,b,None]+Du*g['hu'][:,b,None])
            d2 = vec*(D2*g['h2'][:,b,None]-Du*g['hu'][:,b,None])
            if a == b:
                d1 += P
                d2 += P*signs
            T += (gram(d1,d1,w)+gram(d2,d2,w))/2
            R += gram(d1+d2,d1+d2,w)/2
    return dict(S=S,T=T,V=V,recoil=R,rhs_density=rhs_density,
                terms=terms,signs=signs)


def response(g, m, c, base, rm, cutoff=1e-13):
    z, eta, E = base['zeta'],base['eta'],base['energy']
    X, info = support(rm['S'],cutoff)
    A = X.T @ (z*z*(rm['T']+eta*rm['recoil'])+z*rm['V']-E*rm['S']) @ X
    A = (A+A.T)/2
    # Each response function is expressed with the scaled dimensionless
    # vector. Only physical D on the RHS introduces a 1/zeta factor.
    d = rm['rhs_density'].T @ (g['w']*(m['P']@c))/z
    dx = X.T @ d
    gaps, U = np.linalg.eigh(A)
    if gaps[0] <= 0:
        raise ValueError('nonpositive response gap')
    strengths = (U.T@dx)**2
    sol = X @ np.linalg.solve(A, dx)
    out = dict(polarizability=float(2*sum(strengths/gaps)),
               dipole_strength_sum=float(sum(strengths)),
               energy_weighted_sum=float(strengths@gaps),
               commutator_reference=1+2*eta,
               lowest_dipole_gap=float(gaps[0]),support=info,
               gaps=gaps.tolist(),strengths=strengths.tolist())
    out['C6_identical'] = dispersion(out, out)
    out['imaginary_frequency_response'] = [
        {'omega':om, 'response':float(2*np.sum(gaps*strengths/(gaps*gaps+om*om)))}
        for om in (0,.1,.25,.5,1,2,5,10,20,50,100)]
    return out, sol


def dispersion(a,b):
    """Second-order nonretarded dipole interaction; no overlap or exchange."""
    ea,eb=np.asarray(a['gaps']),np.asarray(b['gaps'])
    sa,sb=np.asarray(a['strengths']),np.asarray(b['strengths'])
    return float(6*np.sum(sa[:,None]*sb[None,:]/(ea[:,None]+eb[None,:])))


def hydrogen_response(degree=10, cutoff=1e-13):
    """Exact ground radial state; complete l=1 multiplets, analytic integrals.

    The two radial polynomials r and r^2 already span the exact static
    first-order correction. Higher orders improve the frequency response.
    """
    I=lambda n: math.factorial(n)/2.**(n+1)
    N=degree+1
    S=np.zeros((N,N)); H=np.zeros_like(S)
    d=np.array([2/math.sqrt(3)*I(i+4) for i in range(N)])
    for i in range(N):
        for j in range(N):
            S[i,j]=I(i+j+4)
            H[i,j]=.5*((i+1)*(j+1)+2)*I(i+j+2)-(i+j+2)*I(i+j+3)/2+.5*I(i+j+4)-I(i+j+3)
    X,info=support(S,cutoff)
    A=X.T@(H+.5*S)@X
    gaps,U=np.linalg.eigh((A+A.T)/2)
    strengths=(U.T@X.T@d)**2
    out=dict(polarizability=float(2*sum(strengths/gaps)),
             dipole_strength_sum=float(sum(strengths)),
             energy_weighted_sum=float(strengths@gaps),commutator_reference=.5,
             gaps=gaps.tolist(),strengths=strengths.tolist(),support=info,
             degree=degree,occupied_only_response=0.)
    out['C6_identical']=dispersion(out,out)
    return out


def radial_profile(base, radii, q_radial=48, q_angle=64):
    """One-entry radial probability, normalized to one, in a_*^-1.

    Independent radial/angle quadrature, not a histogram of energy nodes.
    This is an envelope distribution, not a trajectory or a sharp boundary.
    """
    y,wy=np.polynomial.laguerre.laggauss(q_radial)
    mu,wm=np.polynomial.legendre.leggauss(q_angle)
    r2,cos=np.meshgrid(y/2,mu,indexing='ij')
    weights=wy[:,None]*y[:,None]**2*wm[None,:]
    zeta=base['zeta']; coeff=np.asarray(base['coefficients'])
    out=[]
    for r in radii:
        rho=zeta*r
        s,t=rho+r2,rho-r2
        u=np.sqrt(np.maximum(0,rho*rho+r2*r2-2*rho*r2*cos))
        P=np.zeros_like(s)
        for c,(l,j,n) in zip(coeff,base['terms']):
            P+=c*s**l*t**j*u**n
        out.append(float(zeta*rho*rho*np.exp(-2*rho)*np.sum(weights*P*P)/16))
    return np.asarray(out)


def solve(degree=6,q=20,eta=0,cutoff=1e-13,response_degree=None):
    g = grid(q)
    m = scalar_matrices(g,degree)
    base,c = ground(m,eta,cutoff)
    base.update(degree=degree,quadrature=q,moments=moments(g,m,c,base['zeta']),
                coefficients=c.tolist(),terms=m['terms'])
    if response_degree is not None:
        rm=response_matrices(g,response_degree)
        base['response'],sc=response(g,m,c,base,rm,cutoff)
        base['response'].update(degree=response_degree,coefficients=sc.tolist(),
                                terms=rm['terms'],signs=rm['signs'].tolist())
    return base


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--degree',type=int,default=6)
    ap.add_argument('--quadrature',type=int,default=20)
    ap.add_argument('--eta',type=float,default=0)
    ap.add_argument('--cutoff',type=float,default=1e-13)
    ap.add_argument('--response-degree',type=int)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    start=time.monotonic()
    result=solve(args.degree,args.quadrature,args.eta,args.cutoff,args.response_degree)
    result['runtime_s']=time.monotonic()-start
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    summary={k:v for k,v in result.items() if k not in ('coefficients','terms','response')}
    if 'response' in result:
        summary['response']={k:v for k,v in result['response'].items()
                             if k not in ('coefficients','terms','signs','gaps','strengths')}
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
