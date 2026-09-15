#!/usr/bin/env python3
"""Silicon middle-plate example: geometry, terminals, helicity, and reactions.

Not a POAMS-derived constitutive law or a measurement of the user's apparatus.
An ideal nested planar guide carries a rigid conducting silicon sheet. Two
clamped-guided beams per axis supply its stiffness; their inertia, coating,
fringing, cross-capacitance and guide cross-coupling are neglected. Patterned
stationary electrodes on the outer supports provide four differential segments
per face. The maintained rotor shafts are a separate subsystem, not assumed to
excite this mode. All dimensions, Q and reference capacitance are declared.

SI throughout. C(g,q)=C_ref*g_ref/g*(1+s*q_alpha/overlap) is the local
parallel-overlap approximation, not a full manufactured capacitance solution.
Fixed terminal voltages require the source-inclusive cost -C V^2/2. This is
not the fixed-charge boundary condition for an isolated charged middle plate.
"""
from dataclasses import asdict, dataclass, replace
import json
import math
import numpy as np


@dataclass(frozen=True)
class Spec:
    young: float = 130e9       # Pa, <100> silicon: IUCr 2014 material comparison
    density: float = 2329.1    # kg/m^3, nominal 20 C silicon: NISTIR 6969
    side: float = 1e-3
    thickness: float = 20e-6
    beam_width: float = 4e-6
    length_x: float = 500e-6
    length_y: float = 510e-6
    quality: float = 200.      # chosen loss parameter, NOT measured
    capacitance: float = 10e-15  # per segment at reference gap, assumed calibration target
    reference_gap: float = 10e-6
    overlap: float = 100e-6
    bias: float = 10.          # central sheet at circuit reference; face segments DC biased
    amplitude: float = .1      # peak, not rms, small AC modulation per segment


def mechanical(spec=Spec()):
    m = spec.density*spec.side**2*spec.thickness
    lengths = np.array([spec.length_x, spec.length_y])
    # One guided beam: 12 E I/L^3 with I=t*w^3/12. Two beams per axis.
    k = 2*spec.young*spec.thickness*spec.beam_width**3/lengths**3
    natural = np.sqrt(k/m)
    damping = m*natural/spec.quality
    return m, k, damping, natural


def rotation(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]])


def circular_basis():
    return np.array([[1., 1.], [1j, -1j]])/math.sqrt(2)


def response(frequency=None, phase=math.pi/2, ratio=1., gaps=None, spec=Spec(), angle=0.):
    m, k, damping, natural = mechanical(spec)
    if frequency is None:
        frequency = natural[0]/(2*math.pi)
    if gaps is None:
        gaps = (spec.reference_gap, spec.reference_gap)
    if frequency <= 0 or min(gaps) <= 0 or spec.quality <= 0:
        raise ValueError('Positive frequency, gaps and quality required.')
    gaps = np.asarray(gaps)
    omega = 2*math.pi*frequency
    r = rotation(angle)
    stiffness = r@np.diag(k)@r.T
    loss = r@np.diag(damping)@r.T
    susceptibility = np.linalg.inv(stiffness-m*omega**2*np.eye(2)-1j*omega*loss)
    voltage = spec.amplitude*np.array([[1, 1j], ratio*np.exp(1j*phase)*np.array([1, -1j])])
    c = spec.capacitance*spec.reference_gap/gaps
    gain = 2*c*spec.bias/spec.overlap
    force_phasors = gain[:, None]*voltage
    a = susceptibility@force_phasors.sum(axis=0)
    # Axial attraction magnitudes, top and bottom. Subtract DC analytically
    # before computing the small incremental force to avoid cancellation loss.
    increment = c/gaps*(.5*np.sum(abs(voltage)**2, axis=1)
                       + spec.bias/spec.overlap*np.real(voltage.conj()@a))
    dc = 2*c/gaps*spec.bias**2
    axial = np.array([-increment[0], increment[0]-increment[1], increment[1]])
    full_axial = np.array([-dc[0], dc[0]-dc[1], dc[1]])+axial
    powers = .5*np.real(force_phasors.conj()@(-1j*omega*a))
    dissipated = .5*omega**2*np.vdot(a, loss@a).real
    momentum = m*omega*np.imag(a[0].conjugate()*a[1])
    kinetic = .25*m*omega**2*np.vdot(a, a).real
    elastic = .25*np.vdot(a, stiffness@a).real
    basis = circular_basis()
    circular_response = basis.conj().T@susceptibility@basis
    return dict(frequency=frequency, omega=omega, amplitude=a, voltages=voltage,
                force_phasors=force_phasors, susceptibility=susceptibility,
                circular_response=circular_response, mass=m, stiffness=stiffness,
                damping=loss, gap_attraction_increment=increment, dc_attraction=dc,
                axial_increment=axial, full_axial=full_axial, powers=powers,
                dissipation=dissipated, momentum=momentum, kinetic=kinetic,
                elastic=elastic, gaps=gaps, axis_angle=angle)


def frozen_cycle_cost(a, voltages, positions, spec=Spec()):
    """AC excess of source-inclusive cycle cost, fixed a and terminal voltages."""
    gaps = np.array([positions[0]-positions[1], positions[1]-positions[2]])
    c = spec.capacitance*spec.reference_gap/gaps
    return -np.sum(c*(.5*np.sum(abs(voltages)**2, axis=1)
                      +spec.bias/spec.overlap*np.real(voltages.conj()@a)))


def terminal_cycle(result, spec=Spec(), count=8192):
    """Independent real-time voltage/charge calculation, including C_dot."""
    w = result['omega']
    t = np.arange(count)*2*math.pi/(count*w)
    carrier = np.exp(-1j*w*t)
    q = np.real(carrier[:, None]*result['amplitude'])
    qdot = np.real(-1j*w*carrier[:, None]*result['amplitude'])
    v = np.real(carrier[:, None, None]*result['voltages'])
    vdot = np.real(-1j*w*carrier[:, None, None]*result['voltages'])
    c0 = spec.capacitance*spec.reference_gap/result['gaps']
    terminal_power = np.zeros((count, 2))
    attraction = np.zeros((count, 2))
    lateral = np.zeros((count, 2, 2))
    for p in range(2):
        for alpha in range(2):
            for sign in [1, -1]:
                c = c0[p]*(1+sign*q[:, alpha]/spec.overlap)
                cdot = c0[p]*sign*qdot[:, alpha]/spec.overlap
                u = spec.bias+sign*v[:, p, alpha]
                udot = sign*vdot[:, p, alpha]
                current = cdot*u+c*udot
                terminal_power[:, p] += u*current
                attraction[:, p] += .5*c*u*u/result['gaps'][p]
                lateral[:, p, alpha] += .5*c0[p]*sign/spec.overlap*u*u
    return dict(q=q, qdot=qdot, terminal_power=terminal_power,
                attraction=attraction, lateral=lateral)


def serial(result):
    keys = ['frequency', 'momentum', 'kinetic', 'elastic', 'dissipation']
    out = {k: float(result[k]) for k in keys}
    for key in ['axial_increment', 'full_axial', 'powers', 'dc_attraction']:
        out[key] = result[key].tolist()
    out['amplitude_real'] = result['amplitude'].real.tolist()
    out['amplitude_imag'] = result['amplitude'].imag.tolist()
    return out


def verify():
    checks = []
    def check(label, passed):
        if not passed:
            raise AssertionError(label)
        checks.append(label)
    def near(a, b, scale=1., tol=2e-8):
        return np.max(abs(np.asarray(a)-np.asarray(b))) <= tol*scale

    spec = Spec()
    m, k, gamma, natural = mechanical(spec)
    # Independent beam energy from the guided cubic shape, not a repeated k formula.
    s = np.linspace(0, 1, 10001)
    displacement = 1e-8
    inertia = spec.thickness*spec.beam_width**3/12
    for j, length in enumerate([spec.length_x, spec.length_y]):
        curvature = displacement*(6-12*s)/length**2
        both_beam_energy = spec.young*inertia*np.trapezoid(curvature**2, s*length)
        check(f'axis {j}: bending energy recovers guide stiffness', near(both_beam_energy, .5*k[j]*displacement**2, .5*k[j]*displacement**2, 3e-8))
    beam_fraction = 2*spec.density*spec.thickness*spec.beam_width*(spec.length_x+spec.length_y)/m
    check('neglected beam mass below one percent of ideal sheet', beam_fraction < .01)
    check('example is kHz, not GHz or THz', 1000 < min(natural)/(2*math.pi) < max(natural)/(2*math.pi) < 1500)

    # Rotationally symmetric linear response may split helicities, not mix them.
    basis = circular_basis()
    antisymmetric = np.array([[0., -1.], [1., 0.]])
    symmetric_response = (.3+.2j)*np.eye(2)+(.07-.03j)*antisymmetric
    helicity = basis.conj().T@symmetric_response@basis
    check('isotropic plus gyrotropic response has no helicity conversion', max(abs(helicity[0,1]), abs(helicity[1,0])) < 1e-15)
    check('gyrotropic response can distinguish helicities without mixing', abs(helicity[0,0]-helicity[1,1]) > .01)
    for theta in [2*math.pi/3, math.pi/2, .713]:
        r = rotation(theta)
        check(f'axial symmetry rotation {theta:.3f}: response commutes', near(r@symmetric_response, symmetric_response@r))

    isotropic = replace(spec, length_y=spec.length_x)
    iso0 = response(phase=0, ratio=.7, spec=isotropic)
    for phi in [.5, math.pi/2, 2.7]:
        isot = response(phase=phi, ratio=.7, spec=isotropic)
        check(f'isotropic plate, phase {phi:.2f}: same axial loads', near(isot['axial_increment'], iso0['axial_increment'], 1e-11))
    base = response(spec=spec)
    chi = base['susceptibility']
    check('directional conversion is (chi_x-chi_y)/2', near(base['circular_response'][0,1], (chi[0,0]-chi[1,1])/2, 100.))
    angle = .371
    turned = response(spec=spec, angle=angle)
    check('rotation of directional structure carries twice the angle', near(turned['circular_response'][0,1], base['circular_response'][0,1]*np.exp(-2j*angle), 100.))

    fixtures = []
    cases = [(natural[0]/(2*math.pi), math.pi/2, 1., (1.,1.), spec),
             (1180., .63, .6, (.93,1.07), spec),
             (1168., 2.4, 1.3, (1.08,.92), spec),
             (1200., .8, .8, (1.,1.), isotropic)]
    for j, (f, phase, ratio, gf, sample) in enumerate(cases):
        gaps = np.array(gf)*sample.reference_gap
        result = response(f, phase, ratio, gaps, sample)
        a = result['amplitude']; forces = result['force_phasors']
        check(f'{j}: actual displacement is within linear overlap', max(abs(a))/sample.overlap < .001)
        check(f'{j}: summed axial reactions cancel', near(sum(result['axial_increment']), 0, 1e-11))
        check(f'{j}: mechanical power equals passive dissipation', near(sum(result['powers']), result['dissipation'], 1e-13))
        real = terminal_cycle(result, sample)
        check(f'{j}: full V*d(CV)/dt input gives mechanical power', near(real['terminal_power'].mean(axis=0), result['powers'], 1e-13))
        inc_real = real['attraction'].mean(axis=0)-result['dc_attraction']
        check(f'{j}: real voltage-cycle axial force equals phasor result', near(inc_real, result['gap_attraction_increment'], 1e-11, 5e-8))
        prescribed = np.real(np.exp(-1j*result['omega']*np.arange(8192)*2*math.pi/(8192*result['omega']))[:,None,None]*forces)
        check(f'{j}: segment voltages give specified lateral force', near(real['lateral'], prescribed, 1e-9))
        j_sample = m*np.mean(real['q'][:,0]*real['qdot'][:,1]-real['q'][:,1]*real['qdot'][:,0])
        check(f'{j}: actual material trajectory carries computed AM', near(j_sample, result['momentum'], 1e-20))
        z = np.array([gaps[0], 0., -gaps[1]])
        derivative = []
        step = 1e-10
        for axis in range(3):
            dz = np.eye(3)[axis]*step
            derivative.append(-(frozen_cycle_cost(a,result['voltages'],z+dz,sample)-frozen_cycle_cost(a,result['voltages'],z-dz,sample))/(2*step))
        check(f'{j}: force signs equal frozen source-inclusive derivative', near(derivative, result['axial_increment'], 1e-11))
        check(f'{j}: common translation does not change internal cost', near(frozen_cycle_cost(a,result['voltages'],z+.01,sample), frozen_cycle_cost(a,result['voltages'],z,sample), 1e-16))
        drive_torque = .5*np.real(a[0]*forces[:,1].conj()-a[1]*forces[:,0].conj()).sum()
        kt = -result['stiffness']@a
        gt = 1j*result['omega']*result['damping']@a
        restoring = .5*np.real(a[0]*kt[1].conj()-a[1]*kt[0].conj())
        damping_torque = .5*np.real(a[0]*gt[1].conj()-a[1]*gt[0].conj())
        check(f'{j}: torque includes the anisotropic guide and damping reactions', near(drive_torque+restoring+damping_torque, 0, 1e-17))
        fixtures.append(dict(frequency=f, phase=phase, ratio=ratio, gap_factors=gf,
                             isotropic=sample.length_x==sample.length_y, output=serial(result)))

    for phi in [0., .3, math.pi/2, math.pi, 1.5*math.pi]:
        result = response(phase=phi)
        predicted = -4*spec.capacitance**2*spec.bias**2*spec.amplitude**2/(spec.reference_gap*spec.overlap**2)*np.imag(chi[0,0]-chi[1,1])*math.sin(phi)
        check(f'equal drives phase {phi:.3f}: signed analytic axial result', near(result['axial_increment'][1], predicted, 1e-12))
    reversed_phase = response(phase=1.5*math.pi)
    check('phase reversal reverses center force', near(reversed_phase['axial_increment'][1], -base['axial_increment'][1], 1e-12))
    check('default phase gives negative center force, not lift', base['axial_increment'][1] < 0)
    stopped_bias = response(spec=replace(spec,bias=0))
    check('zero bias means no fundamental lateral actuation in this layout', max(abs(stopped_bias['amplitude'])) == 0)
    check('zero bias equal signals still exert balanced axial pressure', max(abs(stopped_bias['full_axial'])) > 0 and stopped_bias['full_axial'][1] == 0)
    double_signal = response(spec=replace(spec,amplitude=2*spec.amplitude))
    check('signal doubled: displacement doubles and incremental force quadruples', near(double_signal['amplitude'], 2*base['amplitude'], 1e-8) and near(double_signal['axial_increment'], 4*base['axial_increment'], 1e-11))
    check('same-pan source and electrode reactions do not add a DC support force', sum(base['full_axial']) == 0.)
    # At equal axes chi_x=chi_y; with real (lossless off-resonance) chi their
    # difference cannot make this particular antisymmetric axial observable.
    real_chi = 1/(k-m*(2*math.pi*1000.)**2)
    check('equal-drive phase load requires dissipative quadrature in this geometry', np.imag(real_chi[0]-real_chi[1]) == 0.)
    return dict(checks_passed=len(checks), checks=checks, units='SI',
                specimen_status='specified ideal example, not fabricated or calibrated in fact',
                inputs=asdict(spec), stiffness_N_per_m=k.tolist(),
                natural_frequencies_Hz=(natural/(2*math.pi)).tolist(),
                neglected_total_beam_mass_fraction=beam_fraction,
                default=serial(base), opposite_phase=serial(reversed_phase), fixtures=fixtures,
                unresolved=['measured full capacitance matrix and its derivatives',
                            'actual damping, guide/coating inertia and cross-coupling',
                            'rotor shaft to driven-mode coupling',
                            'GHz/THz material realization', 'Earth-relative joint action'])


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
