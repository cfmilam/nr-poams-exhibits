#!/usr/bin/env python3
"""Phase -> local angular momentum -> complete support ledger.

Conditional continuation of vortex_material_transition.py. The rotating-control
response, inductance and gyromagnetic calibration are INPUTS, not POAMS-derived
material constants. The example is not a measured specimen or a weight claim.
Only Python + NumPy are needed; no target-force fit or domain multiplier.
"""
import json
import math

import numpy as np

from vortex_material_transition import fold, well, v, gauss_integral

CHECKS = []
G, C = 6.67430e-11, 299792458.0
HBAR = 1.0545718176461565e-34
# Declared conventional metrology example, not a microscopic derivation.
PHASE_UNIT = 3.291059784e-16  # Wb per rad of the chosen flux-bias chart
GYRO = -8.7941e10            # (A m^2)/(J s), signed orbital calibration
RING_RADIUS, SPACING = 5e-3, 1e-3
L_MINUS = 30e-9             # H, specified opposed-mode effective inductance
AREA = math.pi * RING_RADIUS**2
KM = PHASE_UNIT**2 / L_MINUS
ZETA, B = 0.5, 2.0
K = KM / (1 - ZETA)
MATRIX = K * np.array([[1., ZETA], [ZETA, 1.]])
SIGMA = AREA / (GYRO * PHASE_UNIT)  # s: du/dOmega, independent rotation probe


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def near(a, b, rtol=1e-9, atol=0.0):
    # Do not let an absolute tolerance of 1 hide errors in tiny SI quantities.
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))


def energy(q, u, matrix=MATRIX):
    x = np.asarray(q) - np.asarray(u)
    return float(x @ matrix @ x / 2 + B * KM * np.sum(1 - np.cos(q)))


def current(q, u, matrix=MATRIX):
    return matrix @ (np.asarray(q) - np.asarray(u)) / PHASE_UNIT


def angular(q, u, matrix=MATRIX):
    return AREA * current(q, u, matrix) / GYRO


def drive_force(q, u, slopes, matrix=MATRIX):
    return float((np.asarray(q) - np.asarray(u)) @ matrix @ slopes)


def receiver_force(z_device, z_receiver, js, receiver_j=1.):
    zcores = z_device + np.array([-SPACING/2, SPACING/2])
    r = z_receiver - zcores
    if np.any(r <= 0):
        raise ValueError("This axial chart has the receiver above both rings")
    return 6 * G * receiver_j / C**2 * float(np.sum(js / r**4))


def receiver_energy(z_device, z_receiver, js, receiver_j=1.):
    r = z_receiver - (z_device + np.array([-SPACING/2, SPACING/2]))
    return -2 * G * receiver_j / C**2 * float(np.sum(js / r**3))


def main():
    qc, up, down = fold(B)
    qp = well(up, B, 1)
    qlow = well(down, B, 0)
    qb, qa, uu = np.array([qc, -qc]), np.array([qp, -qp]), np.array([up, -up])
    ib, ia = current(qb, uu), current(qa, uu)
    jb, ja = angular(qb, uu), angular(qa, uu)
    di, dj = ia - ib, ja - jb
    dq = qp - qc
    inductance = PHASE_UNIT**2 * np.linalg.inv(MATRIX)
    check("positive inductance eigenvalues", np.all(np.linalg.eigvalsh(inductance) > 0))
    check("opposed inductance is declared L_minus",
          near(inductance @ np.array([1., -1.]), L_MINUS * np.array([1., -1.])))
    check("ring energy and calibrated inductive energy agree",
          near(.5 * ia @ inductance @ ia, .5 * (qa-uu) @ MATRIX @ (qa-uu)))
    check("phase stationarity and joining current agree", near(ia[0], -B*KM/PHASE_UNIT*math.sin(qp)))
    check("opposed current steps cancel", di[0] == -di[1])
    check("opposed local angular steps cancel", dj[0] == -dj[1])
    check("signed angular-current conversion", near(dj[0], AREA*di[0]/GYRO))
    check("angular step is not a single hbar", abs(dj[0]/HBAR) > 1e10)

    # Conjugate rotations, not dB/dq (which has energy/torque units).
    domega = 1e-5
    jrot = []
    for i in range(2):
        shift = np.eye(2)[i] * SIGMA * domega
        jrot.append(-(energy(qa, uu+shift)-energy(qa, uu-shift))/(2*domega))
    check("independent rotation derivatives give angular momenta", near(jrot, ja, 1e-9))
    check("common rotation sees zero net opposed angular momentum", near(sum(jrot), 0., atol=1e-34))
    # Integer branch/bias translation preserves material current and this
    # reduced energy. It need not preserve the external source's full state.
    shift = np.array([2*math.pi, -4*math.pi])
    check("integer branch/bias translation preserves current", near(current(qa+shift, uu+shift), ia))
    check("integer branch/bias translation preserves reduced energy", near(energy(qa+shift, uu+shift), energy(qa, uu)))
    # Reversing both contour labels changes components, not physical vectors.
    basis = np.diag([-1., 1.])
    transformed = current(basis@qa, basis@uu, basis@MATRIX@basis)
    check("contour relabel preserves physical current vector", near(basis@transformed, ia))

    # First-ring switch maps to equal/opposite mechanical angular impulses in
    # its support/drive ledger. A torque is not a vertical force without geometry.
    local_support_angular_impulse = -dj
    check("each angular impulse has opposite local receipt", np.all(dj+local_support_angular_impulse == 0))

    # Constant receiver phase slope k=1 rad/m is a reference normalization,
    # not a supplied physical apparatus coupling. Whole balance support = -F.
    slopes = np.array([1., -1.])
    fb, fa = drive_force(qb, uu, slopes), drive_force(qa, uu, slopes)
    df = fa-fb
    dz = 1e-5
    numerical = -(energy(qa, uu+slopes*dz)-energy(qa, uu-slopes*dz))/(2*dz)
    check("receiver position derivative gives signed force", near(numerical, fa, 1e-9))
    check("physical force step equals 2 K_minus k Delta q", near(df, 2*KM*dq))
    check("same derivative via current and bias flux", near(fa, PHASE_UNIT*float(ia@slopes)))
    check("no positional coupling gives no drive force", drive_force(qa, uu, np.zeros(2)) == 0)
    # Put BOTH source and ring on the balance; common translation has no effect.
    def attached_energy(z_d, z_c):
        return energy(qa, uu+slopes*(z_d-z_c))
    fc = -(attached_energy(0,dz)-attached_energy(0,-dz))/(2*dz)
    check("external source reaction is opposite", near(fc, -fa, 1e-9))
    check("source on balance cancels common-translation load",
          attached_energy(.12,.12) == attached_energy(0,0))

    # Whole slow cycle: equal-time linear up/down control ramps, not merely
    # subtracting the force jumps. Integrate each distinct stable branch.
    du = up-down
    low_integral = gauss_integral(lambda q: (q-(q+B*math.sin(q)))*(1+B*math.cos(q)), qlow, qc)
    high_integral = gauss_integral(lambda q: (q-(q+B*math.sin(q)))*(1+B*math.cos(q)), 2*math.pi-qc, qp)
    # Integral x dt per ramp of one second: each u sweeps same interval.
    xint_low, xint_high = low_integral/du, high_integral/du
    check("independent ramp force integrals cancel", abs(xint_low+xint_high) < 1e-12)
    for u in np.linspace(down+1e-7, up-1e-7, 31):
        xl = well(float(u), B, 0)-u
        uh = 2*math.pi-u
        xh = well(float(uh), B, 1)-uh
        if abs(xl+xh) > 2e-12:
            raise AssertionError("half-cycle symmetry")
    check("entire quasistatic cycle has half-period sign reversal", True)
    # All continuous branch changes and both jumps return J; two jumps alone
    # would omit the drive's angular budget.
    qrev = 2*math.pi-qc
    j_start = SIGMA*KM*(qlow-down)
    j_reverse = SIGMA*KM*(qrev-down)
    increments = [jb[0]-j_start, dj[0], j_reverse-ja[0], j_start-j_reverse]
    check("ramps plus both slips return local J", abs(sum(increments))/abs(dj[0]) < 1e-14)
    dwell = .5
    dwell_impulse = fa*dwell
    check("asymmetric dwell can leave reciprocal impulse", dwell_impulse != 0)
    check("dwell impulse cancels against counterpart", dwell_impulse + (-fa)*dwell == 0)
    pair_heat = energy(qb, uu)-energy(qa, uu)
    work = gauss_integral(lambda q: -2*KM*(q-(q+B*math.sin(q)))*(1+B*math.cos(q)), qlow,qc)
    work += gauss_integral(lambda q: -2*KM*(q-(q+B*math.sin(q)))*(1+B*math.cos(q)), qp,qrev)
    check("SI cycle work equals two releases", near(work, 2*pair_heat, 1e-10))
    check("stationary dwell adds no ideal phase work", 0.*float(ia.sum()) == 0)

    # External G-kernel remains a separate conditional channel; partial
    # position derivatives hold measured endpoint J fixed as in section 5C.
    rcv = .1
    fgb = receiver_force(0,rcv,jb)
    fga = receiver_force(0,rcv,ja)
    eps = 1e-6
    fd = -(receiver_energy(eps,rcv,ja)-receiver_energy(-eps,rcv,ja))/(2*eps)
    frc = -(receiver_energy(0,rcv+eps,ja)-receiver_energy(0,rcv-eps,ja))/(2*eps)
    check("G-channel force matches same energy derivative", near(fd,fga,1e-8))
    check("G-channel receiver gets opposite reaction", near(frc,-fga,1e-8))
    expected_dfg = 6*G/C**2*dj[0]*((rcv+SPACING/2)**-4-(rcv-SPACING/2)**-4)
    check("G-channel signed step from measured Delta J", near(fga-fgb,expected_dfg))
    far = 10.
    exact = (far+SPACING/2)**-4-(far-SPACING/2)**-4
    check("balanced far receiver retains extra spacing/range suppression",
          near(exact,-4*SPACING/far**5,rtol=2e-8))

    # If the SAME displacement also changes u, J_z is not zero. Verify the
    # complete derivative (and the envelope theorem on a re-equilibrated branch)
    # in scaled units; epsG is a numerical probe, NOT a physical enhanced G.
    kmat = np.array([[1., .5], [.5, 1.]])
    sig, epsG, sep, rz = -1.7, .03, 1., 3.
    control_slopes = np.array([.13, -.13])
    def full_book(q, z):
        u = uu + control_slopes*z
        x = np.asarray(q)-u
        js = sig*kmat@x
        rr = rz-z-np.array([-sep/2,sep/2])
        return float(.5*x@kmat@x + np.sum(1-np.cos(q)) - epsG*np.sum(js/rr**3))
    def full_force(q, z):
        x = np.asarray(q)-(uu+control_slopes*z)
        js = sig*kmat@x
        jz = -sig*kmat@control_slopes
        rr = rz-z-np.array([-sep/2,sep/2])
        return float(x@kmat@control_slopes+epsG*np.sum(jz/rr**3+3*js/rr**4))
    epsz = 1e-5
    numerical_all = -(full_book(qa,epsz)-full_book(qa,-epsz))/(2*epsz)
    check("combined positional derivative includes changing J", near(numerical_all,full_force(qa,0),1e-8))
    def relaxed(z):
        q = qa.copy()
        rr = rz-z-np.array([-sep/2,sep/2])
        for _ in range(15):
            grad = kmat@(q-(uu+control_slopes*z))+np.sin(q)-epsG*(sig*kmat).T@(rr**-3)
            q -= np.linalg.solve(kmat+np.diag(np.cos(q)),grad)
        return q, grad
    qstar, grad = relaxed(0)
    check("full coupled branch actually stationary", np.linalg.norm(grad) < 1e-12)
    qzp,_ = relaxed(epsz)
    qzm,_ = relaxed(-epsz)
    branch_force = -(full_book(qzp,epsz)-full_book(qzm,-epsz))/(2*epsz)
    check("full branch derivative matches constrained partial derivative", near(branch_force,full_force(qstar,0),1e-8))

    # Static material branches cannot determine the rotation response by
    # themselves: changing sigma changes J but not B at Omega=0.
    check("rotation slope is an independent constitutive input",
          energy(qa,uu) == energy(qa,uu+np.zeros(2)*SIGMA*7) and
          not near(7*ja,ja))

    # Finite spatial circulation: planar orbits carry AM but Wr integrand=0.
    angles = 2*math.pi*np.arange(128)/128
    pts = np.column_stack((RING_RADIUS*np.cos(angles),RING_RADIUS*np.sin(angles),np.zeros(128)))
    tangents = np.column_stack((-np.sin(angles),np.cos(angles),np.zeros(128)))
    momenta = (ja[0]/(128*RING_RADIUS))*tangents
    total_p = np.sum(momenta,axis=0)
    total_j = np.sum(np.cross(pts,momenta),axis=0)
    check("finite circulating entries have no net translation", np.linalg.norm(total_p) < abs(ja[0]/RING_RADIUS)*1e-14)
    check("finite circulating entries retain local orbital AM", near(total_j[2],ja[0]))
    wrnum = np.einsum('ijk,ijk->ij',np.cross(tangents[:,None,:],tangents[None,:,:]),pts[:,None,:]-pts[None,:,:])
    check("planar centerline writhe numerator is identically zero", np.all(wrnum == 0))
    frame = np.tile([0.,0.,1.],(128,1))
    check("fixed normal ribbon has no twist variation", np.all(np.diff(frame,axis=0) == 0))
    # About Earth: R x P plus local AM. Purely circulating relative momenta
    # do not change R x P and are not an Earth-orbital increment.
    rearth = np.array([6.371e6,0.,0.])
    change_jearth_orb = np.cross(rearth,total_p)
    check("Earth-orbital change requires translational momentum",
          np.linalg.norm(change_jearth_orb) < np.linalg.norm(rearth)*abs(ja[0]/RING_RADIUS)*1e-14)

    # Natural-radius law is a separate constrained-energy derivative.
    # Work in scaled units to check sign independently without roundoff.
    mu, mass, ell, radius = 3.,2.,4.,5.
    def radial_book(r):
        return -mu*mass/r+ell**2/(2*mass*r**2)
    support = mu*mass/radius**2-ell**2/(mass*radius**3)
    eps = 1e-4
    check("radial support sign from constraint derivative", near((radial_book(radius+eps)-radial_book(radius-eps))/(2*eps),support,1e-8))
    # Same local q-energy and J calibration at r0, different radial derivative:
    # an allowed but unselected Earth-coupling function. This is a logical
    # nonidentifiability witness, not a proposed effect.
    state_fn = lambda q: float(np.sum(1-np.cos(q)))
    for coefficient in (-2.,0.,3.):
        extra_at_r0 = coefficient*(radius-radius)*state_fn(qa)
        check(f"unselected radial term invisible to local switch ({coefficient})", extra_at_r0 == 0)
    check("local switch cannot select a unique radial load", state_fn(qa) != state_fn(qb))

    result = {
        "scope":"conditional calibrated ring example; no measured apparatus or new weight law",
        "inputs":{"radius_m":RING_RADIUS,"spacing_m":SPACING,"L_minus_H":L_MINUS,
                  "phase_unit_Wb":PHASE_UNIT,"signed_gyro":GYRO,"rotation_slope_s":SIGMA,"b":B,"zeta":ZETA},
        "K_minus_J":KM,"q_before":qc,"q_after":qp,"u_up":up,"u_down":down,
        "I_before_A":ib.tolist(),"I_after_A":ia.tolist(),"delta_I_A":di.tolist(),
        "J_before_Js":jb.tolist(),"J_after_Js":ja.tolist(),"delta_J_Js":dj.tolist(),
        "delta_J_over_hbar":float(dj[0]/HBAR),"pair_heat_J":pair_heat,"cycle_work_J":work,
        "force_step_per_bias_gradient_J":df,"force_before_per_bias_gradient_J":fb,
        "force_after_per_bias_gradient_J":fa,
        "symmetric_cycle_impulse_per_gradient_Js":2*KM*(xint_low+xint_high),
        "extra_upper_dwell_s":dwell,"dwell_impulse_per_gradient_Js":dwell_impulse,
        "G_force_step_N_per_receiver_Js_at_point1m":fga-fgb,
        "checks_passed":len(CHECKS),"checks":CHECKS,
    }
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
