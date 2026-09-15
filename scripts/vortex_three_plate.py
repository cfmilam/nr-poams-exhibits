#!/usr/bin/env python3
"""Driven material modes in a three-plate geometry; conditional, not a weight law.

q_i is an actual planar displacement of an effective participating material mass.
This is a small-amplitude constitutive example, NOT a first-principles POAMS
material law, electrical-circuit design, or GHz/THz performance prediction.
All examples are normalized unless explicitly labelled SI. NumPy required.
"""
from dataclasses import dataclass, replace
import json
import math
import numpy as np


@dataclass(frozen=True)
class Model:
    mass: tuple = (1., 1., 1.)
    stiffness: tuple = (1., 1., 1.)
    damping: tuple = (.08, .08, .08)
    coupling: float = .2
    length: float = 1.
    reference_gap: float = 1.
    anisotropy: float = 0.


def matrices(gaps=(1., 1.), model=Model()):
    """At fixed bias, positive springs with kappa(g)=kappa0 exp(-(g-g0)/ell)."""
    if min(gaps) <= 0 or abs(model.anisotropy) >= 1:
        raise ValueError('Require nonzero gaps and |anisotropy|<1.')
    coupling = model.coupling*np.exp(-(np.array(gaps)-model.reference_gap)/model.length)
    kt, kb = coupling
    base = np.diag(model.stiffness)+np.array([[kt, -kt, 0], [-kt, kt+kb, -kb], [0, -kb, kb]])
    delta = np.diag([0., model.stiffness[1]*model.anisotropy, 0.])
    return base+delta, base-delta, coupling


def drives(ratio=1., phase=0., handedness=(1, -1), amplitude=.005):
    """Re[f exp(-i w t)]. Physical helicity is not a scalar voltage phase."""
    f = np.zeros((3, 2), dtype=complex)
    f[0] = amplitude*np.array([1., 1j*handedness[0]])
    f[2] = amplitude*ratio*np.exp(1j*phase)*np.array([1., 1j*handedness[1]])
    return f


def solve(omega=1.1, ratio=1., phase=0., handedness=(1, -1), gaps=(1., 1.), model=Model(), amplitude=.005):
    kx, ky, coupling = matrices(gaps, model)
    f = drives(ratio, phase, handedness, amplitude)
    mass, gamma = np.array(model.mass), np.array(model.damping)
    a = np.empty((3, 2), dtype=complex)
    for axis, k in enumerate([kx, ky]):
        a[:, axis] = np.linalg.solve(k-np.diag(omega**2*mass+1j*omega*gamma), f[:, axis])
    squared = np.sum(abs(a)**2, axis=1)
    differences = a[[0, 1]]-a[[1, 2]]
    gap_load = -.25*coupling/model.length*np.sum(abs(differences)**2, axis=1)
    forces = np.array([-gap_load[0], gap_load[0]-gap_load[1], gap_load[1]])
    kinetic = .25*omega**2*mass*squared
    potential = .25*sum(np.vdot(a[:, j], k@a[:, j]).real for j, k in enumerate([kx, ky]))
    angular_momentum = mass*omega*np.imag(np.conj(a[:, 0])*a[:, 1])
    drive_power = .5*np.real(np.sum(np.conj(f)*(-1j*omega*a), axis=1))
    loss_power = .5*omega**2*gamma*squared
    drive_torque = .5*np.real(a[:, 0]*np.conj(f[:, 1])-a[:, 1]*np.conj(f[:, 0]))
    damping_torque = -gamma*omega*np.imag(np.conj(a[:, 0])*a[:, 1])
    # The central anisotropy is anchored to its plate; the plate takes the opposite torque.
    anisotropy_torque = (kx[1, 1]-ky[1, 1])*.5*np.real(a[1, 0]*np.conj(a[1, 1]))
    return dict(a=a, f=f, forces=forces, gap_load=gap_load, kinetic=kinetic,
                energy=float(np.sum(kinetic)+potential), angular_momentum=angular_momentum,
                drive_power=drive_power, loss_power=loss_power,
                drive_torque=drive_torque, damping_torque=damping_torque,
                anisotropy_torque=float(anisotropy_torque))


def frozen_energy(a, gaps, model=Model()):
    """Cycle average at FIXED amplitudes; do not differentiate a re-solved driven state."""
    kx, ky, _ = matrices(gaps, model)
    return .25*sum(np.vdot(a[:, j], k@a[:, j]).real for j, k in enumerate([kx, ky]))


def gyro_potential(theta, angular_momentum=1., cage_rate=.1, inertia=.5):
    # Routh reduction of the spin-cyclic symmetric rotor, with phi_dot=cage_rate imposed.
    return -.5*inertia*cage_rate**2*np.sin(theta)**2-angular_momentum*cage_rate*np.cos(theta)


def gyro_tilt_torque(theta, angular_momentum=1., cage_rate=.1, inertia=.5):
    return np.sin(theta)*(inertia*cage_rate**2*np.cos(theta)-angular_momentum*cage_rate)


def gyro_trajectory(r=.005, damping=1., theta0=math.pi-.01, step=.002, duration=22.):
    """Dimensionless time sqrt(J*Omega/I_perp)t, with J*Omega>0 and r=I_perp*Omega/J.

    The center is stationary, one tilt is free, azimuth is externally imposed.
    Not a claim that a torque-free three-gimbal mount obeys this constraint.
    """
    def rhs(y):
        th, v = y
        return np.array([v, math.sin(th)*(r*math.cos(th)-1)-damping*v])
    y = np.array([theta0, 0.])
    energy = []
    dissipation = 0.
    crossed = None
    for n in range(round(duration/step)+1):
        energy.append(.5*y[1]**2-.5*r*math.sin(y[0])**2-math.cos(y[0]))
        if crossed is None and y[0] < math.pi/2:
            crossed = n*step
        previous_v = y[1]
        a = rhs(y); b = rhs(y+step*a/2); c = rhs(y+step*b/2); d = rhs(y+step*c)
        y += step*(a+2*b+2*c+d)/6
        dissipation += .5*damping*step*(previous_v**2+y[1]**2)
    final_energy = .5*y[1]**2-.5*r*math.sin(y[0])**2-math.cos(y[0])
    return y, energy, dissipation, final_energy, crossed


def verify():
    checks = []
    def check(label, condition):
        if not condition:
            raise AssertionError(label)
        checks.append(label)
    def near(a, b, atol=2e-10):
        return np.allclose(a, b, rtol=2e-9, atol=atol)

    symmetric = solve()
    imbalanced = solve(ratio=.6)
    check('matched opposed pair gives zero central incremental load', abs(symmetric['forces'][1]) < 1e-14)
    check('unequal amplitudes give nonzero center load', abs(imbalanced['forces'][1]) > 1e-8)
    kx, _, _ = matrices()
    check('three collective resonances: 1,sqrt(1.2),sqrt(1.6)', near(np.linalg.eigvalsh(kx), [1., 1.2, 1.6]))
    check('total internal load cancels for unequal amplitudes', near(sum(imbalanced['forces']), 0))

    for index, (w, ratio, phase, hs, gaps, anis) in enumerate([
        (1.1, 1., 0., (1, -1), (1., 1.), 0.),
        (.7, .4, 1.3, (1, -1), (.8, 1.2), 0.),
        (1.23, .7, .9, (1, 1), (1.2, .8), .3),
        (1.04, .6, -1.2, (-1, 1), (.9, 1.1), .4),
    ]):
        model = replace(Model(), anisotropy=anis)
        result = solve(w, ratio, phase, hs, gaps, model)
        a, f = result['a'], result['f']
        check(f'{index}: power input equals material loss', near(sum(result['drive_power']), sum(result['loss_power'])))
        check(f'{index}: torque includes anisotropy mount', near(sum(result['drive_torque']+result['damping_torque'])+result['anisotropy_torque'], 0))
        check(f'{index}: total three-plate force zero', near(sum(result['forces']), 0))
        kx, ky, _ = matrices(gaps, model)
        check(f'{index}: positive stiffness matrices', min(np.linalg.eigvalsh(kx))>0 and min(np.linalg.eigvalsh(ky))>0)
        z = np.array([gaps[0], 0., -gaps[1]])
        def energy_at(positions):
            return frozen_energy(a, (positions[0]-positions[1], positions[1]-positions[2]), model)
        derivative = []
        for j in range(3):
            dz = np.eye(3)[j]*1e-5
            derivative.append(-(energy_at(z+dz)-energy_at(z-dz))/2e-5)
        check(f'{index}: three force signs from frozen-coordinate derivative', near(derivative, result['forces']))
        check(f'{index}: common translation leaves energy unchanged', near(energy_at(z+19.), energy_at(z)))
        # Integrate real trajectories over a cycle independently of phasor expressions.
        times = np.arange(4096)*2*np.pi/(4096*w)
        q = np.real(a[None, :, :]*np.exp(-1j*w*times[:, None, None]))
        qdot = np.real(-1j*w*a[None, :, :]*np.exp(-1j*w*times[:, None, None]))
        actual_j = np.mean(q[:, :, 0]*qdot[:, :, 1]-q[:, :, 1]*qdot[:, :, 0], axis=0)*np.array(model.mass)
        check(f'{index}: angular momentum matches sampled actual circulation', near(actual_j, result['angular_momentum']))
        differences = q[:, [0, 1], :]-q[:, [1, 2], :]
        couplings = matrices(gaps, model)[2]
        real_gap_load = -.5*couplings/model.length*np.mean(np.sum(differences*differences, axis=2), axis=0)
        check(f'{index}: force average matches sampled cycle', near(real_gap_load, result['gap_load']))

    reference = solve(ratio=.6, gaps=(.8, 1.2))
    for phase in [.31, 1.7, math.pi]:
        changed = solve(ratio=.6, gaps=(.8, 1.2), phase=phase)
        check(f'opposite isotropic circular drives: phase {phase} leaves loads unchanged', near(changed['forces'], reference['forces']))
    anis = replace(Model(), anisotropy=.3)
    a0 = solve(ratio=.6, model=anis, phase=0.)
    a1 = solve(ratio=.6, model=anis, phase=1.2)
    check('specified anisotropy enables phase dependence', not near(a0['forces'], a1['forces'], atol=1e-9))
    reversed_senses = solve(ratio=.6, handedness=(-1, 1))
    check('helicity reversal preserves isotropic gap loads', near(reversed_senses['forces'], imbalanced['forces']))
    check('helicity reversal reverses AM without reversing gap load', near(reversed_senses['angular_momentum'], -imbalanced['angular_momentum']))
    twice = solve(ratio=.6, amplitude=.01)
    check('doubling force amplitude quadruples stored energy, loads and loss', near([twice['energy'], *twice['forces'], sum(twice['loss_power'])], 4*np.array([imbalanced['energy'], *imbalanced['forces'], sum(imbalanced['loss_power'])])))

    # Mass-center response: center plate motion with a held chassis is a transient scale signal.
    phase = np.arange(4096)*2*np.pi/4096
    acceleration = -2.3**2*.01*np.cos(phase)
    check('periodic internal displacement: mean total support increment zero', abs(np.mean(.4*acceleration)) < 1e-15)
    moving_mass, frame_mass = .4, 2.
    frame_shift = -moving_mass/(frame_mass+moving_mass)*.02
    check('free chassis recoil retains complete COM', near(frame_mass*frame_shift+moving_mass*(frame_shift+.02), 0))

    for th in [.1, .8, 2.4]:
        derivative = -(gyro_potential(th+1e-6)-gyro_potential(th-1e-6))/2e-6
        check(f'gyro {th}: imposed-rate Routh potential differentiates correctly', near(derivative, gyro_tilt_torque(th)))
    check('reversing imposed cage rate reverses preferred alignment', near(gyro_potential(.2, cage_rate=.1), gyro_potential(math.pi-.2, cage_rate=-.1)))
    check('no imposed cage rate: no alignment potential in this model', gyro_potential(.8, cage_rate=0.) == 0.)
    y, energies, heat, end, crossed = gyro_trajectory()
    check('with damping imposed-azimuth rotor approaches co-spin', abs(y[0]) < .002)
    check('gyro effective energy decreases', max(np.diff(energies)) < 1e-12)
    check('gyro effective energy decrease equals damping integral', abs(energies[0]-end-heat)<3e-6)
    earth_rate = 7.292115e-5
    cage_rate = 2*math.pi/60  # ONE rpm, illustration only.
    spin_rate = 10000*2*math.pi/60
    j_over_i = 2*spin_rate  # Ideal thin disk: I_axis/I_diameter=2, not a measured gyro.
    alignment_time = 1/math.sqrt(j_over_i*cage_rate-cage_rate*cage_rate)
    check('slow one-rpm cage exceeds terrestrial rotation rate by >1000', cage_rate/earth_rate > 1000)
    check('illustrative characteristic alignment time under .1s', alignment_time < .1)

    output = {
        'checks_passed': len(checks), 'checks': checks,
        'normalization': 'm0=k0=ell=1; omega0=sqrt(k0/m0); forces in k0*ell, NOT newtons',
        'illustrative_imbalanced_case': {key: value.tolist() if isinstance(value, np.ndarray) else value
                                       for key, value in imbalanced.items() if key not in ['a','f']},
        'gyro_illustration_not_fitted': {'cage_rpm': 1, 'spin_rpm':10000, 'cage_to_earth_rate_ratio':cage_rate/earth_rate,
                                        'small_tilt_characteristic_s':alignment_time,
                                        'dimensionless_half_turn_crossing_time':crossed},
        'unselected': ['actual material and terminal-to-mode force mapping', 'bias-dependent stiffness/coupling/loss',
                       'shaft rotation to fast-mode selection', 'Earth-relative joint interaction'],
    }
    return output


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
