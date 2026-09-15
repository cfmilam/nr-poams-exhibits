#!/usr/bin/env python3
"""Finite angular exchange in the carried Vortex Cores pair model.

Run: python3 scripts/vortex_discrete_exchange.py
Requires NumPy. No measured force or target transition rate is fitted.

MODEL INPUTS, not consequences of integer linking:
  * A fixed SU(2) angular multiplet at each specified circulation.
  * B_ij = -(2G/c^2) L_i.L_j / r_ij^3, used as an operator.
  * i*hbar*da/dt = B*a, as in the corpus's finite configuration ledger.
  * Fixed relational positions (a constrained, not freely recoiling, chart).

The labels m below are angular projections, NOT geometric Lk, Tw or Wr.
Evolution gives coherent configuration weights, NOT a stochastic reconnection
or reset law. The three-circulation example uses j=1, not a macroscopic device.
Its endpoint forces are the carried interaction's contributions; holding forces
have the opposite signs. Common rest-energy constants are omitted once.
"""

import itertools
import json
import math

import numpy as np

H = 6.62607015e-34
HBAR = H / (2 * math.pi)
C = 299792458.0
G = 6.67430e-11  # one calibrated coupling, CODATA 2022 central value


def assert_near(actual, expected, atol=2e-11):
    assert np.allclose(actual, expected, rtol=2e-11, atol=atol), (actual, expected)


def angular_matrices(j):
    """Dimensionless L/hbar, ordered m=j,j-1,...,-j."""
    assert j >= 0 and float(2 * j).is_integer()
    m = np.arange(j, -j - 1, -1, dtype=float)
    plus = np.zeros((len(m), len(m)), complex)
    for col, value in enumerate(m):
        if col:
            plus[col - 1, col] = math.sqrt(j * (j + 1) - value * (value + 1))
    minus = plus.conj().T
    return m, [(plus + minus) / 2, (plus - minus) / (2j), np.diag(m)]


def embed(matrix, address, sizes):
    result = np.ones((1, 1), complex)
    for index, size in enumerate(sizes):
        result = np.kron(result, matrix if index == address else np.eye(size))
    return result


def angular_book(js):
    local = [angular_matrices(j) for j in js]
    sizes = [len(item[0]) for item in local]
    basis = list(itertools.product(*(item[0] for item in local)))
    operators = [[embed(matrix, i, sizes) for matrix in item[1]] for i, item in enumerate(local)]
    pairs = {(i, k): sum(operators[i][a] @ operators[k][a] for a in range(3))
             for i, k in itertools.combinations(range(len(js)), 2)}
    return basis, operators, pairs


def sector(basis, operators, total_m):
    indices = [i for i, entry in enumerate(basis) if abs(sum(entry) - total_m) < 1e-12]
    return [basis[i] for i in indices], {key: matrix[np.ix_(indices, indices)]
                                      for key, matrix in operators.items()}


def pair_book(j):
    basis, _, pairs = angular_book([j, j])
    basis, pairs = sector(basis, pairs, 0)
    return basis, -pairs[0, 1]


def pair_tridiagonal(j):
    m = np.arange(j, -j - 1, -1)
    matrix = np.diag(m * m)
    for index, value in enumerate(m[:-1]):
        entry = -(j * (j + 1) - value * (value - 1)) / 2
        matrix[index, index + 1] = matrix[index + 1, index] = entry
    return matrix


def eigensystem(matrix, initial):
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    coefficients = eigenvectors.conj().T @ initial
    return eigenvalues, eigenvectors, coefficients


def state_at(system, time):
    values, vectors, coefficients = system
    return vectors @ (np.exp(-1j * values * time) * coefficients)


def expected(state, operator):
    return float(np.vdot(state, operator @ state).real)


def integrated_expectation(system, operator, time):
    """Exact spectral integral from 0 to time, including degenerate levels."""
    values, vectors, coefficients = system
    gaps = values[:, None] - values[None, :]
    phase_integral = time * np.exp(0.5j * gaps * time) * np.sinc(gaps * time / (2 * np.pi))
    amplitudes = coefficients.conj()[:, None] * coefficients[None, :]
    return float(np.sum(amplitudes * (vectors.conj().T @ operator @ vectors) * phase_integral).real)


def spatial_book(positions, pairs, disabled=()):
    """Positions in units d; B in E0=2G*hbar^2/(c^2*d^3); F in E0/d."""
    template = next(iter(pairs.values()))
    bill = np.zeros_like(template)
    forces = [np.zeros_like(template) for _ in positions]
    for (i, k), matrix in pairs.items():
        if (i, k) in disabled:
            continue
        delta = positions[i] - positions[k]
        bill -= matrix / abs(delta) ** 3
        force = -3 * matrix * delta / abs(delta) ** 5
        forces[i] += force
        forces[k] -= force
    return bill, forces


def main():
    report = {"scope": __doc__, "checks": []}

    # Independent tensor-product vs ladder/tridiagonal construction.
    for j in (0.5, 1, 2, 3):
        m, (x, y, z) = angular_matrices(j)
        assert_near(x @ y - y @ x, 1j * z)
        assert_near(x @ x + y @ y + z @ z, j * (j + 1) * np.eye(len(m)))
        basis, bill = pair_book(j)
        assert_near(bill, pair_tridiagonal(j))
        total_j = np.arange(0, 2 * j + 1)
        spectrum = -0.5 * (total_j * (total_j + 1) - 2 * j * (j + 1))
        assert_near(np.linalg.eigvalsh(bill), sorted(spectrum))
        # The pair book is a polynomial in the total Casimir.
        total_squared = 2 * j * (j + 1) * np.eye(len(m)) - 2 * bill
        assert_near(bill @ total_squared - total_squared @ bill, np.zeros_like(bill))
        initial = np.zeros(len(m))
        initial[0] = 1
        fluctuation = (bill - j * j * np.eye(len(m))) @ initial
        assert_near(np.linalg.norm(fluctuation), j)
        system = eigensystem(bill, initial)
        for time in (0.001, 0.01, 0.1, 0.3, 1):
            state = state_at(system, time)
            outside_weight = float(sum(abs(state[1:]) ** 2))
            assert outside_weight <= min(1, (j * time) ** 2) + 1e-12
        for value in m[:-1]:
            delta_omega_in_units_kappa_hbar = -(2 * value - 1)
            cost = lambda v: v * v + delta_omega_in_units_kappa_hbar * v
            assert_near(cost(value), cost(value - 1))
    report["checks"].append("SU(2), Casimir spectrum, four full pair blocks, diagonal crossing formula, and finite-time departure bound")

    # An exact j=1 opposed pair: no two-state truncation and no irreversible jump.
    basis, bill = pair_book(1)
    assert_near(bill, [[1, -1, 0], [-1, 0, -1], [0, -1, 1]])
    initial = np.array([1.0, 0, 0])
    system = eigensystem(bill, initial)
    pair_states = []
    for time in (0, 0.1, 1, math.pi / 3, 2, 10):
        state = state_at(system, time)
        assert_near(np.vdot(state, state), 1)
        assert_near(expected(state, bill), 1)
        assert_near(abs(state[1]) ** 2, 4 * math.sin(1.5 * time) ** 2 / 9)
        # Radial force is 3B/d; its expectation is unchanged by B evolution.
        assert_near(expected(state, 3 * bill), 3)
        pair_states.append({"tau": time, "weights": (abs(state) ** 2).tolist(), "B_over_E0": expected(state, bill)})
    report["pair_j1"] = {"basis": basis, "B_over_E0": bill.real.tolist(), "snapshots": pair_states,
                         "central_max_weight": 4 / 9, "first_central_max_tau": math.pi / 3,
                         "emitted_energy_in_isolated_model_J": 0}
    report["checks"].append("exact j=1 analytic central weight, conserved energy, unchanged internal constraint")

    # A fully specified external receiver C: j=1 at all three addresses.
    full_basis, ops, full_pairs = angular_book([1, 1, 1])
    positions = np.array([-0.5, 0.5, 3.0])
    full_bill, _ = spatial_book(positions, full_pairs)
    for axis in range(3):
        total = sum(item[axis] for item in ops)
        assert_near(full_bill @ total - total @ full_bill, np.zeros_like(full_bill))
    basis, pairs = sector(full_basis, full_pairs, 0)
    assert len(basis) == 7
    bill, forces = spatial_book(positions, pairs)
    device_force = forces[0] + forces[1]
    assert_near(device_force + forces[2], np.zeros_like(bill))
    assert_near(device_force, 3 * (pairs[0, 2] / 3.5 ** 4 + pairs[1, 2] / 2.5 ** 4))
    # Differentiate the energy independently under common translation of A,B.
    step = 1e-5
    shift = np.array([step, step, 0])
    plus, _ = spatial_book(positions + shift, pairs)
    minus, _ = spatial_book(positions - shift, pairs)
    assert_near(device_force, -(plus - minus) / (2 * step), atol=3e-10)
    reflected_bill, reflected_forces = spatial_book(-positions, pairs)
    assert_near(reflected_bill, bill)
    assert_near(reflected_forces[0] + reflected_forces[1], -device_force)
    decoupled, isolated_forces = spatial_book(positions, pairs, disabled=((0, 2), (1, 2)))
    assert_near(isolated_forces[0] + isolated_forces[1], np.zeros_like(bill))
    report["checks"].append("total vector AM; energy-gradient endpoint forces; reflection sign; receiver-disconnection null")

    initial_index = basis.index((1, -1, 0))
    initial = np.zeros(len(basis))
    initial[initial_index] = 1
    transitions = []
    for target, inverse_cube in [((0, 0, 0), 1), ((0, -1, 1), 3.5 ** -3), ((1, 0, -1), 2.5 ** -3)]:
        index = basis.index(target)
        assert_near(bill[index, initial_index], -inverse_cube)
        transitions.append({"from": [1, -1, 0], "to": target, "matrix_entry_over_E0": -inverse_cube,
                            "diagonal_cost_difference_over_E0": float((bill[index, index] - bill[initial_index, initial_index]).real)})
    system = eigensystem(bill, initial)
    snapshots = []
    for time in (0, 0.1, 1, 10, 20):
        state = state_at(system, time)
        weights = abs(state) ** 2
        assert_near(sum(weights), 1)
        assert_near(expected(state, bill), 1)
        assert_near(sum(sum(entry) * weight for entry, weight in zip(basis, weights)), 0)
        assert_near(expected(state, device_force), -expected(state, forces[2]))
        internal_bill = -pairs[0, 1]
        external_bill = bill - internal_bill
        snapshots.append({"tau": time,
                          "receiver_nonzero_m_weight": float(sum(weight for entry, weight in zip(basis, weights) if entry[2] != 0)),
                          "receiver_mean_Lz_over_hbar": float(sum(entry[2] * weight for entry, weight in zip(basis, weights))),
                          "B_AB_over_E0": expected(state, internal_bill),
                          "B_AC_plus_BC_over_E0": expected(state, external_bill),
                          "B_total_over_E0": expected(state, bill),
                          "device_interaction_force_over_E0_per_d": expected(state, device_force),
                          "receiver_interaction_force_over_E0_per_d": expected(state, forces[2]),
                          "integrated_device_interaction_force_over_hbar_per_d": integrated_expectation(system, device_force, time)})
    # Independent quadrature of the analytic spectral integral (dimensionless).
    times = np.linspace(0, 1, 4001)
    sampled = [expected(state_at(system, time), device_force) for time in times]
    quadrature = np.trapezoid(sampled, times) if hasattr(np, "trapezoid") else np.trapz(sampled, times)
    assert_near(integrated_expectation(system, device_force, 1), quadrature, atol=3e-10)
    # A unitary chart change must not change the same physical bill or response.
    phases = np.exp(1j * np.arange(len(basis)) * 0.37)
    U = np.diag(phases)
    transformed_system = eigensystem(U @ bill @ U.conj().T, U @ initial)
    transformed_force = U @ device_force @ U.conj().T
    assert_near(expected(state_at(transformed_system, 1), transformed_force), snapshots[2]["device_interaction_force_over_E0_per_d"])
    report["checks"].append("specified exchange entries; norm/energy conservation; reciprocal response; spectral integral; chart invariance")
    report["receiver_example"] = {"j_each": 1, "positions_over_d": positions.tolist(), "initial": [1, -1, 0],
                                  "basis": basis, "transitions": transitions,
                                  "B_over_E0": bill.real.tolist(), "snapshots": snapshots,
                                  "interpretation": "Coherent exchange weights and fixed-position interaction contributions, not isolated calorimetric bursts or stochastic event times."}

    # The physical units remain those of the weak carried kernel, not pump frequency.
    d = 0.1  # declared illustration: not a physical j=1 apparatus specification
    kappa = 2 * G / (C ** 2 * d ** 3)
    E0 = kappa * HBAR ** 2
    t0 = 1 / (kappa * HBAR)
    F0 = E0 / d
    report["illustrative_SI_conversion"] = {"d_m": d, "E0_J": E0, "t0_s": t0, "F0_N": F0,
                                          "tau1_device_force_N": snapshots[2]["device_interaction_force_over_E0_per_d"] * F0,
                                          "tau1_receiver_force_N": snapshots[2]["receiver_interaction_force_over_E0_per_d"] * F0,
                                          "tau1_integrated_interaction_Ns": snapshots[2]["integrated_device_interaction_force_over_hbar_per_d"] * HBAR / d}
    report["initial_elementary_mixing_scales_not_event_rates"] = [
        {"assumed_max_projection_J_s": J, "g_s_inverse": kappa * J, "inverse_g_s": 1 / (kappa * J),
         "departure_weight_upper_bound_at_one_second": min(1, (kappa * J) ** 2)}
        for J in (1e3, 1e5, 1e7)]
    # Restoring a prescribed state against its evolved correlations is not a free
    # reset. This script introduces no detector, dissipative bath or pump reset.
    report["result"] = "PASS: explicit conditional exchange and counterpart calculated; no topological trigger, domain burst or laboratory sensitivity claim"
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
