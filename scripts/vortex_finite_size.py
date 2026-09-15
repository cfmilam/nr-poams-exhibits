#!/usr/bin/env python3
"""Finite-size continuations of the carried Earth-relative mixed-AM model.

This does NOT select a microscopic POAMS interaction. Two explicit finite sums
recover the same point-body law and survive co-located subdivision, yet give
different separated-pair responses. A third, point-invisible tensor term shows
why the horizontal null is not universal. All are first-order, fixed-rate books.

Positions refer to Earth; J_i refer to named closed subcirculations, not an
arbitrary partition into open arcs. All participating mount/drive changes must
be included. The apparatus boundary is fixed, not an adjustable normalization.
The kernel's derivatives are derivatives of a relational chart, not a medium.
Requires Python + NumPy and the neighboring material/phase/Earth scripts.
"""
from decimal import Decimal, localcontext
import json
import math

import numpy as np

from vortex_earth_constraint import A, OMEGA_E, R_E
from vortex_material_transition import fold, well, gauss_integral
from vortex_phase_support import G, C, KM, SIGMA, B

CHECKS = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def near(x, y, rtol=1e-9, atol=1e-12):
    return bool(np.allclose(x, y, rtol=rtol, atol=atol))


def derivative(fn, value=0.):
    return np.imag(fn(value + 1e-25j)) / 1e-25


def bvec(x, axis):
    """b = [axis - (axis.u)u]/r, analytic also for complex-step checks."""
    rr = np.sqrt(np.dot(x, x))
    return np.cross(x, np.cross(axis, x)) / rr**3


def dbvec(x, axis, direction):
    rr = np.sqrt(np.dot(x, x))
    xu, zu, zx = np.dot(x, direction), np.dot(axis, direction), np.dot(axis, x)
    return -(axis*xu + x*zu + direction*zx)/rr**3 + 3*zx*x*xu/rr**5


def cost(kind, x, angular, y, masses, axis, coefficient=1.):
    """coefficient = a_G*Omega; C1 and C2 are alternative, not additive, laws."""
    total_mass, total_j = np.sum(masses), np.sum(angular, axis=0)
    if kind == 'circulation':
        return -coefficient*total_mass*sum(np.dot(j, bvec(xx, axis)) for xx, j in zip(x, angular))
    if kind == 'orbital':
        return -coefficient*np.dot(total_j, sum(m*bvec(yy, axis) for yy, m in zip(y, masses)))
    raise ValueError(kind)


def support(kind, x, angular, y, masses, axis, direction, coefficient=1.):
    """Common rigid displacement, fixed J, bias/calibration and orientation."""
    total_mass, total_j = np.sum(masses), np.sum(angular, axis=0)
    if kind == 'circulation':
        return -coefficient*total_mass*sum(np.dot(j, dbvec(xx, axis, direction)) for xx, j in zip(x, angular))
    if kind == 'orbital':
        return -coefficient*np.dot(total_j, sum(m*dbvec(yy, axis, direction) for yy, m in zip(y, masses)))
    raise ValueError(kind)


def first_moment(x, angular, center):
    return sum(np.outer(xx-center, j) for xx, j in zip(x, angular))


def moment_witness(x, angular, y, masses, axis, coefficient=1., chi=1.):
    """Not adopted: allowed finite tensor term, invisible to point-body data."""
    mass = np.sum(masses)
    center = np.sum(masses[:, None]*y, axis=0)/mass
    rr = np.sqrt(np.dot(center, center))
    q = first_moment(x, angular, center)
    return chi*coefficient*mass*np.dot(axis, center)*np.trace(q)/rr**3


def coaxial(radius, spacing, tilt, latitude, amplitude=1.):
    # Orthonormal frame: e0=radial, e1=north; geographic latitude fixes axis.
    radial = np.array([1., 0., 0.])
    axis = np.array([math.sin(latitude), math.cos(latitude), 0.])
    normal = np.array([math.cos(tilt), math.sin(tilt), 0.])
    center = radius*radial
    x = np.array([center+spacing*normal/2, center-spacing*normal/2])
    angular = np.array([amplitude*normal, -amplitude*normal])
    return x, angular, axis, radial, normal


def gamma_factor(normal, radial, axis):
    s, c, t = np.dot(axis, radial), np.dot(normal, radial), np.dot(axis, normal)
    return s*(3*c*c-1)-2*c*t


def decimal_tilted_load(amplitude, spacing='0.001'):
    """Exact two-center C1 derivative at equator/45deg, using 65 digits.

    Avoid subtracting two ~1e-68 N terms at ordinary precision to quote a
    ~1e-77 N difference. Material amplitude precision is still illustrative.
    """
    with localcontext() as ctx:
        ctx.prec = 65
        r, d = Decimal(str(R_E)), Decimal(spacing)
        j = Decimal(str(amplitude))
        unit = Decimal(1)/Decimal(2).sqrt()
        a = 2*Decimal(str(G))/Decimal(str(C))**2
        k = a*Decimal(str(OMEGA_E))

        def projected_derivative(sign):
            x, y = r+sign*d*unit/2, sign*d*unit/2
            norm = (x*x+y*y).sqrt()
            return -unit*(x+y)/norm**3 + 3*x*y*unit*(x+y)/norm**5

        exact = -k*j*(projected_derivative(1)-projected_derivative(-1))
        leading = -2*k*j*d/r**3
        return {'exact_N': str(exact), 'leading_N': str(leading),
                'relative_finite_correction': str((exact-leading)/leading)}


def main():
    rng = np.random.default_rng(81273)
    axis = np.array([0., 0., 1.])
    x = np.array([[3.1, .2, .7], [2.9, -.1, .8]])
    y = np.array([[3., .1, .6], [3.2, -.2, .7]])
    angular = np.array([[.7, -.2, .4], [-.1, .3, .2]])
    masses = np.array([.4, .8])
    direction = np.array([.8, .3, -.1])
    check('b equals the point-body Xi/r vector',
          near(bvec(x[0], axis), (axis-np.dot(axis, x[0])*x[0]/np.dot(x[0], x[0]))/np.linalg.norm(x[0])))
    check('analytic kernel gradient matches complex-step vector derivative',
          near(dbvec(x[0], axis, direction), derivative(lambda t: bvec(x[0]+t*direction, axis))))

    # All candidate costs are scalar fixed-rate books, with explicit derivatives.
    for kind in ['circulation', 'orbital']:
        fn = lambda t: cost(kind, x+t*direction, angular, y+t*direction, masses, axis)
        n = support(kind, x, angular, y, masses, axis, direction)
        check(kind+' support matches independent cost differentiation', near(n, derivative(fn)))
        check(kind+' Earth reaction is equal/opposite', near(-n, derivative(lambda t: cost(kind, x-t*direction, angular, y-t*direction, masses, axis))))
        check(kind+' common translation leaves the relative book unchanged',
              near(cost(kind, (x+direction)-direction, angular, (y+direction)-direction, masses, axis), fn(0.)))
        # Any first-order cost -Omega*B has a compatible fixed-rate Hamiltonian.
        omega, inertia = .37, 2.8
        b0 = -cost(kind, x, angular, y, masses, axis)
        def h(p, coupling): return p*p/(2*inertia)-p*coupling*b0/inertia
        def legendre(coupling):
            p = inertia*omega+coupling*b0
            return h(p, coupling)-omega*p
        check(kind+' source-work Legendre completion preserves first-order cost',
              near(derivative(legendre), -omega*b0))
        # Physical locations fixed: no partial sum m_i*J_i normalization.
        xs = np.concatenate([x[:1], x[:1], x[1:]])
        js = np.concatenate([angular[:1]*.3, angular[:1]*.7, angular[1:]])
        ys = np.concatenate([y[:1], y[:1], y[1:]])
        ms = np.array([masses[0]*.3, masses[0]*.7, masses[1]])
        check(kind+' invariant under co-located circulation subdivision', near(cost(kind, xs, js, y, masses, axis), fn(0.)))
        check(kind+' invariant under co-located mass subdivision', near(cost(kind, x, angular, ys, ms, axis), fn(0.)))
        check(kind+' labels/permutation cannot change support', near(support(kind,x[::-1],angular[::-1],y[::-1],masses[::-1],axis,direction), n))
        rot = np.linalg.qr(rng.normal(size=(3,3)))[0]
        check(kind+' invariant under a common proper rotation', near(cost(kind,x@rot.T,angular@rot.T,y@rot.T,masses,rot@axis),fn(0.)))
        common = np.repeat([[3., .2, .7]],2,axis=0)
        point = -np.sum(masses)*np.dot(np.sum(angular,axis=0),bvec(common[0],axis))
        check(kind+' recovers the unsplit point-body law',near(cost(kind,common,angular,common,masses,axis),point))

    # Matched physical pair, rigidly separated along the common ring normal.
    radius, spacing, tilt, latitude = 3., .015, .6, .4
    xx,jj,zz,radial,normal = coaxial(radius,spacing,tilt,latitude)
    mass_pair = np.array([.5,.5])
    n1 = support('circulation',xx,jj,xx,mass_pair,zz,radial)
    n2 = support('orbital',xx,jj,xx,mass_pair,zz,radial)
    check('two coarsening-consistent extensions give different finite loads', abs(n1)>1e-6 and n2==0.)
    q = first_moment(xx,jj,radius*radial)
    check('first spatial AM moment is J*d*n tensor n',near(q,spacing*np.outer(normal,normal)))
    check('matched first moment is independent of chosen expansion origin',near(q,first_moment(xx,jj,radius*radial+np.array([.2,.4,-.1]))))
    gam = gamma_factor(normal,radial,zz)
    check('coaxial factor equals the meridian tilt formula',near(gam,-math.sin(latitude)*math.sin(tilt)**2-math.cos(latitude)*math.sin(2*tilt)))
    leading = 2*spacing*gam/radius**3
    check('exact finite derivative agrees with leading spatial moment', near(n1,leading,rtol=1e-4))
    xhalf,jhalf,_,_,_ = coaxial(radius,spacing/2,tilt,latitude)
    nhalf = support('circulation',xhalf,jhalf,xhalf,mass_pair,zz,radial)
    check('next matched-pair load term is cubic in spacing',near((n1-leading)/(nhalf-leading/2),8.,rtol=2e-4))
    # Geometry factor is bounded by 2/sqrt(3), independently of latitude.
    factors=[]
    for lat in np.linspace(-math.pi/2,math.pi/2,41):
        s,k=math.sin(lat),math.cos(lat)
        eig=np.linalg.eigvalsh(np.array([[0.,-k,0.],[-k,-s,0.],[0.,0.,-s]]))
        factors.extend(abs(eig))
    check('coaxial geometry has an order-one bound, not an amplifier',max(factors)<=2/math.sqrt(3)+1e-12)
    s=1/math.sqrt(3)
    check('the geometric bound is sharp',near((s+math.sqrt(4-3*s*s))/2,2/math.sqrt(3)))
    for degrees in [-70,0,30,80,90]:
        xh,jh,zh,rh,_=coaxial(3.,.2,0.,math.radians(degrees))
        check('horizontal coaxial C1 radial load vanishes at latitude '+str(degrees),
              near(support('circulation',xh,jh,xh,mass_pair,zh,rh),0.))
    check('C2 matched null holds at arbitrary separation/orientation', n2==0.)
    check('reversing both circulation senses reverses the finite C1 load',near(support('circulation',xx,-jj,xx,mass_pair,zz,radial),-n1))

    # An independent permitted tensor term prevents overclaiming a universal null.
    xh,jh,zh,rh,_=coaxial(3.,.2,0.,math.pi/6)
    def witness_shift(t): return moment_witness(xh+t*rh,jh,xh+t*rh,mass_pair,zh)
    wn=derivative(witness_shift)
    check('point-invisible tensor term changes horizontal support',near(wn,-2*.2*math.sin(math.pi/6)/3**3) and wn!=0.)
    check('witness is linear in the explicitly unselected coefficient',near(moment_witness(xh,jh,xh,mass_pair,zh,chi=-1),-witness_shift(0.)))
    xc,jc,zc,_,_=coaxial(3.,0.,0.,math.pi/6)
    check('witness disappears in the co-located point limit',moment_witness(xc,jc,xc,mass_pair,zc)==0.)
    xref=np.concatenate([xh[:1],xh[:1],xh[1:]])
    jref=np.concatenate([jh[:1]/2,jh[:1]/2,jh[1:]])
    yref=np.concatenate([xh[:1],xh[:1],xh[1:]])
    check('witness also survives co-located subdivision',near(moment_witness(xref,jref,yref,np.array([.25,.25,.5]),zh),witness_shift(0.)))
    rot=np.linalg.qr(rng.normal(size=(3,3)))[0]
    check('witness has scalar rotational symmetry',near(moment_witness(xh@rot.T,jh@rot.T,xh@rot.T,mass_pair,rot@zh),witness_shift(0.)))
    check('witness has spatial parity symmetry (J and Earth axis are axial)',near(moment_witness(-xh,jh,-xh,mass_pair,zh),witness_shift(0.)))
    check('witness has time-reversal symmetry (Omega and J both reverse)',near(moment_witness(xh,-jh,xh,mass_pair,zh,coefficient=-1),witness_shift(0.)))

    # Whole ramps as well as jumps: force is linear in actual J(t), not event rate.
    qc,up,down=fold(B)
    qp=well(up,B,1)
    dj=SIGMA*KM*(qp-qc)
    ramp1=gauss_integral(lambda u:well(u,B,0)-u,down,up,180)
    ramp2=gauss_integral(lambda u:well(u,B,1)-u,down,up,180)
    check('equal-speed full switch/reset ramps cancel the linear finite response',abs(ramp1+ramp2)<1e-12)
    check('a fixed off-axis dwell can leave reciprocal impulse',n1*(qp-up)!=0.)
    check('horizontal C1 geometry stays null for any scalar J history',near(support('circulation',xh,7*jh,xh,mass_pair,zh,rh),0.))

    dec=decimal_tilted_load(dj)
    check('high-precision exact SI residual has the derived positive sign',Decimal(dec['exact_N'])>0)
    check('millimetre SI residual agrees with the moment result',abs(Decimal(dec['relative_finite_correction']))<Decimal('1e-18'))
    check('calibrated AM step is inherited without a new population factor',near(dj,-3.22649888477e-23,rtol=1e-10,atol=0.))
    print(json.dumps({
        'scope':'conditional finite-size continuations, not selected POAMS physics',
        'checks_passed':len(CHECKS),'checks':CHECKS,
        'assumed_total_mass_kg':1.,'assumed_core_spacing_m':.001,
        'assumed_Earth_radius_m':R_E,'nominal_rotation_rad_s':OMEGA_E,
        'calibrated_ring_A_step_Js':dj,
        'C1_horizontal_coaxial_step_N':0.,'C2_any_matched_pair_step_N':0.,
        'C1_equator_45deg_coaxial_step':dec,
        'spacing_over_radius':.001/R_E,
        'tilted_residual_to_uncompensated_north_reference':2*.001/R_E,
        'unselected_horizontal_witness':'-2 chi a_G M Omega sin(latitude) J d/r^3; chi NOT fitted or adopted',
        'remaining':'actual joint spatial coupling and mapping of participating closed circulations; no apparatus-independent no-go'
    },indent=2))


if __name__=='__main__':
    main()
