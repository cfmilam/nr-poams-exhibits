#!/usr/bin/env python3
"""Earth-relative support from the carried mixed-AM Hamiltonian.

Conditional continuation, not a first-principles derivation of that interaction.
The chosen body-level coupling is H_mix = -a J_star.L/r^3, a=2G/c^2.
J_star's identification with a material's local circulation is an explicit
additional assumption. Fixed Earth rotation is enforced by H-Omega*p; p is
canonical orbital action, not automatically mechanical m*r^2*Omega.

All geometric/derivative tests use dimensionless inputs to expose differences.
SI results use declared illustrative inputs, never a target-force fit.
Only Python + NumPy and the two neighboring material/phase scripts are used.
"""
import json
import math

import numpy as np

from vortex_material_transition import fold, well, gauss_integral
from vortex_phase_support import G, C, KM, SIGMA, B

CHECKS = []
A = 2 * G / C**2
R_E = 6.371e6              # illustrative spherical radius, m
MU_E = 3.986004418e14       # calibrated reference GM, m^3/s^2
OMEGA_E = 7.292115e-5       # nominal mean rotation, rad/s


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def near(x, y, rtol=1e-10, atol=0.):
    return bool(np.allclose(x, y, rtol=rtol, atol=atol))


def deriv(fn, x):
    """Complex-step derivative of analytic test functions, no subtractive loss."""
    return np.imag(fn(x + 1e-25j)) / 1e-25


def hamiltonian(r, p, jn, mass, mu, latitude, a):
    s = math.cos(latitude)
    inertia = mass * r**2 * s**2
    return p**2/(2*inertia) - mu*mass/r - a*jn*p/(s*r**3)


def p_stationary(r, omega, jn, mass, latitude, a):
    s = math.cos(latitude)
    return mass*r**2*s**2*(omega + a*jn/(s*r**3))


def fixed_rate_book(r, omega, jn, mass, mu, latitude, a):
    # Exact for the stipulated truncated Hamiltonian. Its a^2 terms are NOT
    # promoted to second-order physical predictions of an uncompleted model.
    s = math.cos(latitude)
    return -mu*mass/r - .5*mass*r**2*s**2*(omega+a*jn/(s*r**3))**2


def fixed_rate_support(r, omega, jn, mass, mu, latitude, a):
    s = math.cos(latitude)
    return (mu*mass/r**2-mass*r*omega**2*s**2
            +a*mass*omega*s*jn/r**2+2*mass*a*a*jn*jn/r**5)


def xi(j, rhat, axis):
    """Projection on Earth's rotation axis perpendicular to the radial ray."""
    return np.dot(j, axis)-np.dot(axis, rhat)*np.dot(j, rhat)


def main():
    # Canonical completion: choose a nonpolar latitude and scaled coefficients.
    r, p, jn, mass, mu, lat, a, omega = 3., 1.7, .6, 1.2, 4., .4, .025, .2
    s = math.cos(lat)
    ps = p_stationary(r, omega, jn, mass, lat, a)
    check("fixed-rate stationarity dH/dp=Omega",
          near(deriv(lambda pp: hamiltonian(r, pp, jn, mass, mu, lat, a), ps), omega))
    check("H-Omega*p equals the completed square",
          near(hamiltonian(r, ps, jn, mass, mu, lat, a)-omega*ps,
               fixed_rate_book(r, omega, jn, mass, mu, lat, a)))
    n = fixed_rate_support(r, omega, jn, mass, mu, lat, a)
    check("full fixed-rate radial derivative",
          near(deriv(lambda rr: fixed_rate_book(rr, omega, jn, mass, mu, lat, a), r), n))
    check("envelope theorem uses fixed p before evaluating p_star",
          near(deriv(lambda rr: hamiltonian(rr, ps, jn, mass, mu, lat, a), r), n))
    check("fixed-action derivative includes coefficient three",
          near(deriv(lambda rr: hamiltonian(rr, p, jn, mass, mu, lat, a), r),
               mu*mass/r**2-p*p/(mass*s*s*r**3)+3*a*jn*p/(s*r**4)))
    check("first-order fixed-rate coefficient is one, not three",
          near(deriv(lambda aa: fixed_rate_support(r, omega, jn, mass, mu, lat, aa), 0.),
               mass*omega*s*jn/r**2))
    check("canonical action differs from mechanical orbital action",
          near(ps-mass*r*r*s*s*omega, a*mass*s*jn/r))
    check("discarded higher order is quadratic in assumed coefficient",
          near(n-(mu*mass/r**2-mass*r*omega**2*s*s+a*mass*omega*s*jn/r**2),
               2*mass*a*a*jn*jn/r**5, rtol=1e-8))
    check("zero coupling recovers held spherical-Earth support",
          near(fixed_rate_support(r, omega, jn, mass, mu, lat, 0.),
               mu*mass/r**2-mass*r*omega**2*s*s))
    # Torque-free circular branch: never differentiate a family of circles as
    # though it were one fixed-action constraint. Positive J shifts r inward.
    l, mass0, mu0, a0, jj = 2., 1., 5., .01, .2
    r0 = l*l/(mass0*mass0*mu0)
    stable_r = (l*l/mass0+math.sqrt((l*l/mass0)**2-12*mu0*mass0*a0*jj*l))/(2*mu0*mass0)
    check("fixed-action circular root has zero constraint",
          abs(deriv(lambda rr: hamiltonian(rr,l,jj,mass0,mu0,0.,a0),stable_r)) < 1e-12)
    check("co-sense circular branch moves inward in the carried model", stable_r < r0)
    def outer_root(j):
        return (l*l/mass0+np.sqrt((l*l/mass0)**2-12*mu0*mass0*a0*j*l))/(2*mu0*mass0)
    check("circular-radius linear shift is -3*a*m*J/L",
          near(deriv(outer_root,0.),-3*a0*mass0/l))

    # Geometry: Earth's axis is not the latitude orbit's full AM direction.
    axis=np.array([0.,0.,1.])
    for degrees in [0.,30.,-47.,89.9]:
        latitude=math.radians(degrees)
        rh=np.array([math.cos(latitude),0.,math.sin(latitude)])
        north=np.array([-math.sin(latitude),0.,math.cos(latitude)])
        position=3.2*rh
        velocity=.13*np.cross(axis,position)
        ell=2.1*np.cross(position,velocity)
        check(f"radial local AM has zero mixed dot product at latitude {degrees}",
              abs(np.dot(rh,ell)) < 1e-14 and abs(xi(rh,rh,axis))<1e-14)
        check(f"north projection carries cos(latitude) at {degrees}",
              near(np.dot(north,ell),2.1*.13*3.2**2*math.cos(latitude)))
        probe=np.array([.8,-.3,.9])
        check(f"coordinate-free projection matches orbital vector at {degrees}",
              near(np.dot(probe,ell),2.1*.13*3.2**2*xi(probe,rh,axis)))
    rh=np.array([1.,0.,0.]); probe=np.array([.8,-.3,.9])
    rot=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
    check("common coordinate rotation leaves Xi unchanged",
          near(xi(rot@probe,rot@rh,rot@axis),xi(probe,rh,axis)))
    check("polar limit vanishes without choosing an orbit plane", xi(probe,axis,axis)==0.)
    check("simultaneous sense reversal preserves dot product",
          near(xi(-probe,rh,-axis),xi(probe,rh,axis)))

    # Matched rings: point-body coupling depends on total participating AM.
    qc, up, down=fold(B)
    qp=well(up,B,1)
    dj=SIGMA*KM*(qp-qc)
    check("calibrated endpoint matches previous phase pass",near(dj,-3.226498885e-23,rtol=1e-9))
    qgrid=np.linspace(qc,qp,41)
    check("opposed participating ring sum cancels through the modeled path",
          all(SIGMA*KM*(q-up)+(-SIGMA*KM*(q-up))==0 for q in qgrid))
    check("housing/drive AM cannot be excluded from body J_star",
          sum([dj,-dj,0.])==0 and sum([dj,0.,-dj])==0)
    # Deliberately unbalanced residual, 1 kg mass: not the matched assembly.
    coeff=A*OMEGA_E/R_E**2
    single=coeff*dj
    check("SI coefficient has positive co-sense support sign", coeff>0 and single<0)
    check("orientation perpendicular to the orbital AM gives zero SI term", coeff*0.==0.)

    # Held controls: Xi(r) need not be constant. Isolate the first-order book
    # to avoid subtracting the enormous baseline from a microscopic correction.
    x0,xslope=.4,-.12
    def x_of_r(rr): return x0+xslope*(rr-r)
    analytic=a*mass*omega*(x0/r**2-xslope/r)
    check("positional calibration contributes the Xi derivative",
          near(deriv(lambda rr: -a*mass*omega*x_of_r(rr)/rr,r),analytic))
    def mutual(zd,ze):
        rr=zd-ze
        return -a*mass*omega*x_of_r(rr)/rr
    nd=deriv(lambda z:mutual(z,0.),r)
    ne=deriv(lambda z:mutual(r,z),0.)
    check("Earth counterpart has opposite mixed positional derivative",near(nd,-ne))
    check("common translation cannot produce a mixed whole-system force",
          near(mutual(r+7.,7.),mutual(r,0.)))

    # Explicit stationary material branch, including its back reaction. Scaled
    # linear Xi=sigma*k*(q-u) suffices; no claimed physical enhancement.
    k0, sig0, uu, eb=.7,.2,1.,.25
    def phase_energy(q, rr):
        return .5*k0*(q-uu)**2+eb*(1-np.cos(q))-a*mass*omega*sig0*k0*(q-uu)/rr
    def relaxed(rr):
        q=.8+0j if isinstance(rr,complex) else .8
        for _ in range(30):
            q-=(k0*(q-uu)+eb*np.sin(q)-a*mass*omega*sig0*k0/rr)/(k0+eb*np.cos(q))
        return q
    qstar=relaxed(r)
    check("full material branch is stationary", abs(deriv(lambda q:phase_energy(q,r),qstar))<1e-14)
    check("relaxed-branch radial derivative equals explicit partial derivative",
          near(deriv(lambda rr:phase_energy(relaxed(rr),rr),r),
               deriv(lambda rr:phase_energy(qstar,rr),r)))

    # Reflection-symmetric hysteresis cycle: this linear residual channel has
    # zero up+down impulse; a paid unequal dwell does not inherit that symmetry.
    qlo=well(down,B,0)
    span=up-down
    def integrand(q): return (q-(q+B*math.sin(q)))*(1+B*math.cos(q))
    ramp_low=gauss_integral(integrand,qlo,qc)/span
    ramp_high=gauss_integral(integrand,2*math.pi-qc,qp)/span
    check("equal-time full ramps cancel the linear AM signal",abs(ramp_low+ramp_high)<1e-12)
    check("phase/reset endpoints return local AM",near(qlo-down,-(qp-up)))
    dwell=.5
    impulse_single=coeff*SIGMA*KM*(qp-up)*dwell
    check("unbalanced residual dwell may have nonzero reciprocal impulse",impulse_single!=0.)
    check("matched rings cancel even time-asymmetric dwell",impulse_single-impulse_single==0.)

    # Topology and local calibration do not select an additional radial law.
    beta=.11
    def witness(rr,q): return beta*(rr-r)*(1-np.cos(q))
    check("unselected radial term leaves local phase energy unchanged",witness(r,qstar)==0.)
    check("unselected radial term leaves local phase gradient unchanged",
          deriv(lambda q:witness(r,q),qstar)==0.)
    check("unselected radial term changes support",near(deriv(lambda rr:witness(rr,qstar),r),beta*(1-np.cos(qstar))))
    # Finite-size extension is not specified by a point-body self-mixed kernel:
    # splitting one body into arbitrary pieces changes sum m_i*J_i.
    check("self-mixed point kernel does not define its own subdivision rule",
          not near(mass*jn,2*(mass/2)*(jn/2)))

    result={
        "scope":"conditional fixed-rate derivative of carried point-body mixed Hamiltonian; no new Earth interaction derived",
        "checks_passed":len(CHECKS),"checks":CHECKS,
        "assumed_radius_m":R_E,"assumed_mu_m3_s2":MU_E,"nominal_omega_rad_s":OMEGA_E,
        "device_mass_kg_for_residual_example":1.,"coupling_a_2G_c2":A,
        "equatorial_load_coefficient_N_per_Js":coeff,
        "calibrated_ring_step_Js":dj,
        "hypothetical_uncompensated_north_step_N":single,
        "matched_point_pair_step_N":0.,"horizontal_ring_mixed_step_N":0.,
        "hypothetical_uncompensated_dwell_impulse_Ns":impulse_single,
        "uncompensated_1e7_Js_equatorial_load_N":coeff*1e7,
        "fixed_action_to_fixed_rate_ratio_same_kinematic_reference":3.,
        "finite_size":"not supplied by the point-body mixed ansatz; no unique multipole coefficient claimed"
    }
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
