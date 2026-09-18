#!/usr/bin/env python3
"""alpha and alpha_G as calibrated dimensionless pair-coupling coefficients.

Companion to the measured-alpha application edition. Uses the same standard-
library Decimal arithmetic and the same conventions as coulomb_calibrated.py.
Computes:
  1. K_Q = alpha * hbar * c                     (EM pair coefficient)
  2. K_G(m_a, m_b) = G * m_a * m_b               (gravitational pair coefficient)
  3. alpha_G(m_a, m_b) = K_G / (hbar * c) = G * m_a * m_b / (hbar * c)
                                                 (dimensionless gravitational pair)
  4. Ratios alpha / alpha_G and the corresponding "gravitational hydrogenic"
     radii and binding-energy scales for a hypothetical purely-gravitational
     hydrogenic system with reduced mass mu = m_a*m_b/(m_a+m_b).

No fitting; no derived origin claim; no cosmological reinterpretation. The
outputs identify two calibrated coefficients and their arithmetical
consequences under the same reciprocal-pair circular model used for the
electromagnetic case in the Coulomb exhibit.
"""
from decimal import Decimal as D, localcontext, getcontext
from pathlib import Path
import hashlib
import json
import platform

HERE = Path(__file__).resolve().parent
SOURCE = "https://physics.nist.gov/cuu/Constants/Table/allascii.txt"

# Inputs: identical to coulomb_calibrated.py so the two runs remain aligned.
INPUTS = {
    "alpha": {"value": "0.0072973525643", "standard_uncertainty": "0.0000000000011",
              "unit": "1", "status": "2022 CODATA measured/adjusted central value"},
    "c": {"value": "299792458", "standard_uncertainty": "0",
          "unit": "m s^-1", "status": "exact SI"},
    "h": {"value": "6.62607015e-34", "standard_uncertainty": "0",
          "unit": "J s", "status": "exact SI"},
    "m_e": {"value": "9.1093837139e-31", "standard_uncertainty": "2.8e-40",
            "unit": "kg", "status": "measured mass calibration; light entry"},
    "m_p": {"value": "1.67262192595e-27", "standard_uncertainty": "5.2e-37",
            "unit": "kg", "status": "measured mass calibration; heavy entry"},
    "G": {"value": "6.67430e-11", "standard_uncertainty": "1.5e-15",
          "unit": "m^3 kg^-1 s^-2",
          "status": "measured input, used only in gravitational branch"},
    "e_SI": {"value": "1.602176634e-19", "standard_uncertainty": "0",
             "unit": "C; also J per eV",
             "status": "exact SI translation only; not an ontology premise"},
}


def two_pi():
    # Chudnovsky pi to current precision, then multiply by 2.
    m, ell, x, k = 1, 13591409, 1, 6
    total = D(ell)
    for j in range(1, getcontext().prec // 14 + 3):
        m = m * (k**3 - 16*k) // j**3
        ell += 545140134
        x *= -262537412640768000
        total += D(m*ell) / D(x)
        k += 12
    pi = 426880 * D(10005).sqrt() / total
    return 2 * pi


def hbar(h):
    return h / two_pi()


def alpha_G_pair(G, m_a, m_b, hbar_val, c):
    """alpha_G(m_a, m_b) = G * m_a * m_b / (hbar * c). Dimensionless."""
    return G * m_a * m_b / (hbar_val * c)


def reduced_mass(m_a, m_b):
    return m_a * m_b / (m_a + m_b)


def hydrogenic_scales(K_pair, mu, hbar_val, c, alpha_dimless):
    """For a reciprocal-pair circular hydrogenic model with pair coefficient K
    (units J·m) and reduced mass mu, and a corresponding dimensionless
    fine-structure ratio alpha_dimless = K / (hbar*c):
        r_1 = hbar / (mu * c * alpha_dimless) = hbar^2 / (mu * K)
        E_1 = -mu * c^2 * alpha_dimless^2 / 2 = -mu * K^2 / (2 * hbar^2)

    K must be given here as its value in J*m. mu and c in SI. alpha_dimless is
    the dimensionless ratio for the pair. Returns (r_1_m, E_1_J).

    For the electromagnetic case, K = alpha*hbar*c and alpha_dimless = alpha.
    For the gravitational case, K = G*m_a*m_b and alpha_dimless = alpha_G.
    """
    r_1_m = hbar_val * hbar_val / (mu * K_pair)
    E_1_J = -mu * K_pair * K_pair / (2 * hbar_val * hbar_val)
    # Cross-check via the alpha form:
    r_1_alt = hbar_val / (mu * c * alpha_dimless)
    E_1_alt = -mu * c * c * alpha_dimless * alpha_dimless / 2
    # These should be equal to double-precision round-off at Decimal precision.
    return {
        "r_1_m": r_1_m, "r_1_m_check": r_1_alt,
        "E_1_J": E_1_J, "E_1_J_check": E_1_alt,
    }


def sha256_source():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def to_str(x):
    return f"{x:.80E}" if isinstance(x, D) else str(x)


def main():
    getcontext().prec = 200
    alpha = D(INPUTS["alpha"]["value"])
    c = D(INPUTS["c"]["value"])
    h = D(INPUTS["h"]["value"])
    m_e = D(INPUTS["m_e"]["value"])
    m_p = D(INPUTS["m_p"]["value"])
    G = D(INPUTS["G"]["value"])
    e_SI = D(INPUTS["e_SI"]["value"])
    hbar_val = hbar(h)

    # Electromagnetic pair coefficient.
    K_Q = alpha * hbar_val * c   # J*m
    # Cross-check dimensionless: alpha_EM = K_Q / (hbar*c) = alpha (by construction).

    pairs = {
        "electron_electron": (m_e, m_e),
        "electron_proton":   (m_e, m_p),
        "proton_proton":     (m_p, m_p),
    }
    results = {}
    for name, (m_a, m_b) in pairs.items():
        a_G = alpha_G_pair(G, m_a, m_b, hbar_val, c)
        K_G = G * m_a * m_b                    # J*m
        mu = reduced_mass(m_a, m_b)

        # Electromagnetic Bohr-analog for this same reduced mass (for ratio):
        em_scales = hydrogenic_scales(K_Q, mu, hbar_val, c, alpha)
        grav_scales = hydrogenic_scales(K_G, mu, hbar_val, c, a_G)

        results[name] = {
            "reduced_mass_kg": str(mu),
            "alpha_G": str(a_G),
            "K_G_J_m": str(K_G),
            "alpha_over_alpha_G": str(alpha / a_G),
            "log10_alpha_over_alpha_G": str(
                (alpha / a_G).ln() / D(10).ln()
            ),
            "em_hydrogenic_r_1_m": str(em_scales["r_1_m"]),
            "em_hydrogenic_E_1_J": str(em_scales["E_1_J"]),
            "em_hydrogenic_E_1_eV": str(em_scales["E_1_J"] / e_SI),
            "grav_hydrogenic_r_1_m": str(grav_scales["r_1_m"]),
            "grav_hydrogenic_E_1_J": str(grav_scales["E_1_J"]),
            "grav_hydrogenic_E_1_eV": str(grav_scales["E_1_J"] / e_SI),
            "r_1_ratio_grav_over_em": str(
                grav_scales["r_1_m"] / em_scales["r_1_m"]
            ),
            "binding_ratio_em_over_grav": str(
                em_scales["E_1_J"] / grav_scales["E_1_J"]
            ),
        }

    # Sanity: for electron_proton, r_1_em should reproduce coulomb_calibrated's
    # finite two-ended Bohr radius 5.2946540946e-11 m within double precision.
    ep = results["electron_proton"]
    r_em_ep_str = ep["em_hydrogenic_r_1_m"]
    r_em_ep = D(r_em_ep_str)
    published_r1 = D("5.29465409461289e-11")
    r_diff_relative = abs(r_em_ep - published_r1) / published_r1

    payload = {
        "schema": "poams.alpha_gravitational_companion.v1",
        "computed_at": "2026-09-17",
        "source_sha256": sha256_source(),
        "source_reference": SOURCE,
        "input_policy": ("No new derivation of alpha or alpha_G. Both are "
                         "calibrated dimensionless pair-coupling coefficients "
                         "in identical structural roles. Uses only measured "
                         "or exact-SI inputs. No cosmological or origin claim."),
        "precision_context": {
            "decimal_prec_significant_digits": 200,
            "reporting_significant_digits": 80,
            "python_version": platform.python_version(),
            "note": ("Input precision is NOT 200 digits; the arithmetic "
                     "precision exceeds the physical precision of the "
                     "inputs. See standard_uncertainty in each input."),
        },
        "inputs": INPUTS,
        "elementary_couplings": {
            "K_Q_J_m": str(K_Q),
            "K_Q_alpha_hbar_c_definition": "alpha * hbar * c",
            "alpha_dimensionless": str(alpha),
        },
        "gravitational_pair_couplings": results,
        "internal_cross_check_electron_proton_r1_vs_coulomb_exhibit": {
            "coulomb_calibrated_finite_two_ended_r_1_m_ref": str(published_r1),
            "this_run_em_r_1_m": r_em_ep_str,
            "relative_agreement": str(r_diff_relative),
            "expected_scope": "within ~1e-11 relative; both derived from the same alpha, m_e, m_p, h, c; any drift signals a broken input link, not a physical result",
        },
        "scope": {
            "what_is_established": [
                "alpha_G(m_a, m_b) is a well-defined dimensionless calibrated pair-coupling coefficient in the same structural role as alpha for the electromagnetic pair.",
                "For every mass pair the analogous reciprocal-pair circular model produces a gravitational-hydrogenic radius and binding-energy scale.",
                "The scale ratio between electromagnetic and gravitational hydrogenic radii equals (alpha/alpha_G).",
            ],
            "what_is_NOT_established": [
                "Numerical value or physical origin of alpha or alpha_G (both accepted as measured or as G*m_a*m_b/hbar*c on measured G and masses).",
                "Existence of a gravitationally-bound hydrogenic atom (the calculated grav radii are typically larger than the observable universe; this is a scale statement, not a physical binding claim).",
                "Any prediction of cosmological structure, dark-matter, or graviton physics.",
                "Any quantitative link between alpha and alpha_G beyond the two calibrations sitting side by side.",
            ],
        },
    }

    output_path = HERE / "alpha_gravitational_companion_results.json"
    output_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {output_path}")

    # Short readable summary to stdout
    print("\n=== Elementary pair coefficients ===")
    print(f"  alpha (measured):     {alpha}")
    print(f"  K_Q = alpha*hbar*c:  {K_Q:.10E} J*m")
    for name, r in results.items():
        print(f"\n=== {name} ===")
        print(f"  alpha_G:                    {D(r['alpha_G']):.10E}")
        print(f"  K_G:                        {D(r['K_G_J_m']):.10E} J*m")
        print(f"  alpha / alpha_G:            {D(r['alpha_over_alpha_G']):.6E}")
        print(f"  log10(alpha/alpha_G):       {float(D(r['log10_alpha_over_alpha_G'])):.4f}")
        print(f"  EM r_1  (finite pair):       {D(r['em_hydrogenic_r_1_m']):.6E} m")
        print(f"  Grav r_1 (finite pair):      {D(r['grav_hydrogenic_r_1_m']):.6E} m")
        print(f"  EM E_1  (finite pair):       {D(r['em_hydrogenic_E_1_eV']):.6E} eV")
        print(f"  Grav E_1 (finite pair):      {D(r['grav_hydrogenic_E_1_eV']):.6E} eV")


if __name__ == "__main__":
    main()
