#!/usr/bin/env python3
"""Measured-alpha back-substitution through existing equations; no fitting.

Standard library only. Decimal arithmetic at 150 and 230 significant digits,
with 80-significant-digit archival outputs. Input precision is NOT 80 digits.
The models below are distinct conditional constructions, not one joint theory.
"""
from decimal import Decimal as D, localcontext, getcontext
from pathlib import Path
import hashlib
import json
import platform

HERE = Path(__file__).resolve().parent
SOURCE = "https://physics.nist.gov/cuu/Constants/Table/allascii.txt"
INPUTS = {
    "alpha": {"value": "0.0072973525643", "standard_uncertainty": "0.0000000000011", "unit": "1", "status": "2022 CODATA measured/adjusted central value"},
    "c": {"value": "299792458", "standard_uncertainty": "0", "unit": "m s^-1", "status": "exact SI"},
    "h": {"value": "6.62607015e-34", "standard_uncertainty": "0", "unit": "J s", "status": "exact SI"},
    "e_SI": {"value": "1.602176634e-19", "standard_uncertainty": "0", "unit": "C; also J per eV", "status": "exact SI translation only; not an elementary ontology premise"},
    "m_light": {"value": "9.1093837139e-31", "standard_uncertainty": "2.8e-40", "unit": "kg", "status": "measured mass calibration; conventional electron-mass correspondence"},
    "m_heavy": {"value": "1.67262192595e-27", "standard_uncertainty": "5.2e-37", "unit": "kg", "status": "measured mass calibration; conventional proton-mass correspondence"},
    "G": {"value": "6.67430e-11", "standard_uncertainty": "1.5e-15", "unit": "m^3 kg^-1 s^-2", "status": "measured input, used only in explicitly gravitational/Planck-conditioned branches"},
}
REFERENCES = {
    "inverse_alpha": "137.035999177",
    "Bohr_radius_m": "5.29177210544e-11",
    "Rydberg_per_m": "10973731.568157",
    "Hartree_J": "4.3597447222060e-18",
}


def atan_series(x):
    """Convergent alternating series, used only for |x| <= 1/5."""
    assert abs(x) <= D("0.2")
    total = term = x
    k = 1
    while True:
        term *= -x*x
        nxt = total + term / (2*k+1)
        if nxt == total:
            return +nxt
        total = nxt
        k += 1


def pi_chudnovsky():
    m, ell, x, k = 1, 13591409, 1, 6
    total = D(ell)
    for j in range(1, getcontext().prec//14 + 3):
        m = m*(k**3-16*k)//j**3
        ell += 545140134
        x *= -262537412640768000
        total += D(m*ell)/D(x)
        k += 12
    return 426880*D(10005).sqrt()/total


def sine(x):
    assert abs(x) < 1
    total = term = x
    j = 1
    while True:
        term *= -x*x/D((2*j)*(2*j+1))
        nxt = total+term
        if nxt == total:
            return +nxt
        total = nxt
        j += 1


def cosine(x):
    assert abs(x) < 1
    total = term = D(1)
    j = 1
    while True:
        term *= -x*x/D((2*j-1)*(2*j))
        nxt = total+term
        if nxt == total:
            return +nxt
        total = nxt
        j += 1


def calculate(precision):
    with localcontext() as ctx:
        ctx.prec = precision
        vals = {k: D(v["value"]) for k,v in INPUTS.items()}
        a, c, h, ev, ml, mh, G = [vals[k] for k in ("alpha", "c", "h", "e_SI", "m_light", "m_heavy", "G")]
        pi = pi_chudnovsky()
        hb = h/(2*pi)
        K = a*hb*c
        mu = ml*mh/(ml+mh)
        rows, checks = {}, []

        def put(key, value, unit, formula, status="algebraic consequence of supplied inputs"):
            rows[key] = {"number": +D(value), "unit": unit, "formula": formula, "status": status}
            return D(value)

        def close(name, x, y, tol=D("1e-110")):
            scale = max(abs(x),abs(y),D("1e-500"))
            error = abs(x-y)/scale
            assert error <= tol, (name, str(error))
            checks.append({"name": name, "passed": True})

        def check(name, condition):
            assert condition, name
            checks.append({"name": name, "passed": True})

        # Separate algorithms for pi; neither uses a stored decimal pi value.
        close("Chudnovsky pi versus independent Machin/arctangent pi", pi, 16*atan_series(D(1)/5)-4*atan_series(D(1)/239))
        put("pi", pi, "1", "Chudnovsky series; independently checked with Machin identity", "mathematical constant, not a fitted parameter")
        put("hbar", hb, "J s", "h/(2*pi)", "fixed SI action unit")
        put("inverse_alpha", 1/a, "1", "1/alpha")
        put("sqrt_alpha", a.sqrt(), "1", "sqrt(alpha)")
        put("inverse_sqrt_alpha", 1/a.sqrt(), "1", "1/sqrt(alpha)")
        put("alpha_squared", a*a, "1", "alpha^2")
        put("binding_fraction", a*a/2, "1", "alpha^2/2", "quadratic circular model: binding/(mu*c^2)")
        put("orbital_to_Compton_period", 1/(a*a), "1", "1/alpha^2", "n=1 circular period divided by h/(mu*c^2); not an event count")
        put("legacy_N", 2/(a*a), "1", "2/alpha^2", "legacy ratio, not a winding integer")
        put("network_kappa", 4*pi*a, "1", "4*pi*alpha", "required weight in the stated network cost; not independently selected")
        put("inverse_network_kappa", 1/(4*pi*a), "1", "1/(4*pi*alpha)", "normalization comparison only")
        put("K_Q", K, "J m", "alpha*hbar*c")
        put("pair_energy_at_1_angstrom", K/D("1e-10")/ev, "eV", "K_Q/(1e-10 m)", "magnitude in the leading 1/r pair law; reference separation only")
        put("pair_force_at_1_angstrom", K/D("1e-10")**2, "N", "K_Q/(1e-10 m)^2", "magnitude in the leading pair law; reference separation only")
        put("relative_circular_rate_v", a*c, "m s^-1", "alpha*c", "n=1 quadratic relative-motion chart, not motion of light")

        # Four-sense energy contrast: independent nuisance energies deliberately added.
        signs = ((1,1),(1,-1),(-1,1),(-1,-1))
        r = D("1e-10")
        ee = [ev*(3+2*s-5*t)+s*t*K/r for s,t in signs]
        dc = sum(D(s*t)*E for (s,t),E in zip(signs,ee))/4
        close("four-sense contrast cancels one-end energies", dc, K/r)
        close("phase acquired over the reference interval r/c", dc/hb*(r/c), a)
        put("connected_phase_per_r_over_c", a, "rad", "(D/hbar)*(r/c)=alpha", "reference exposure, not a derived completion duration or propagation time")
        put("connected_phase_degrees", a*180/pi, "degree", "alpha*180/pi", "phase coordinate, distinct from asin(alpha) cone angle")
        put("full_connected_phase_reference_intervals", 2*pi/a, "1", "2*pi/alpha", "stationary phase-rate comparison only; not an integer count of exchanges")
        put("reference_readout_dark_probability", sine(a/2)**2, "1", "sin^2(alpha/2)", "ideal closed-history readout for c*(I1-I0)=1, st=+1, analysis phase=0; chosen exposure")
        close("closed-history dark and bright outputs sum to unity", sine(a/2)**2+cosine(a/2)**2,D(1))
        close("network back-substitution restores pair coefficient", rows["network_kappa"]["number"]/(4*pi)*hb*c, K)

        # Legacy cone identities, evaluated with cancellation-free small differences.
        b = (1-a*a).sqrt()
        theta = 2*atan_series(a/(1+b))
        cap = 2*pi*a*a/(1+b)
        put("cone_angle", theta, "rad", "asin(alpha)", "legacy cone identification, not a microscopic selection rule")
        put("cone_angle_degrees", theta*180/pi, "degree", "asin(alpha)*180/pi", "legacy cone identification")
        put("cone_dilation_minus_one", a*a/(b*(1+b)), "1", "1/sqrt(1-alpha^2)-1")
        put("cone_cap_solid_angle", cap, "sr", "2*pi*(1-sqrt(1-alpha^2))")
        put("cone_cap_relative_leading_error", a*a/(1+b)**2, "1", "cap/(pi*alpha^2)-1", "nonzero correction to the old leading-order equality")
        close("inverse sine reproduces alpha", sine(theta), a)
        close("exact cap versus stable expression", cap, 2*pi*(1-b))
        close("cap residual formula", cap/(pi*a*a)-1, rows["cone_cap_relative_leading_error"]["number"])

        # Atomic readouts with independently supplied mass calibration.
        put("reduced_mass", mu, "kg", "m_light*m_heavy/(m_light+m_heavy)", "mass inputs; both ends participate")
        for label, mass in (("light_limit", ml), ("two_end", mu)):
            length = hb/(mass*c)
            radius = length/a
            binding = mass*c*c*a*a/2
            R = mass*c*a*a/(2*h)
            omega = mass*c*c*a*a/hb
            period = 2*pi/omega
            put(label+"_Compton_length",length,"m","hbar/(mass*c)","mass-calibrated reference length")
            put(label+"_radius",radius,"m","hbar/(mass*c*alpha)","n=1 quadratic circular model; light_limit uses mass=m_light, two_end uses mass=mu")
            put(label+"_binding_J",binding,"J","mass*c^2*alpha^2/2","gross binding; no fine, radiative or finite-structure corrections")
            put(label+"_binding_eV",binding/ev,"eV","mass*c^2*alpha^2/(2 eV_to_J)","gross binding only")
            put(label+"_binding_mass_equivalent",binding/(c*c),"kg","mass*alpha^2/2","mass-equivalent of the gross spectral limit; not a new constituent mass")
            put(label+"_Rydberg",R,"m^-1","mass*c*alpha^2/(2*h)","gross spectral coefficient; not an independent calibration input")
            put(label+"_limit_frequency",binding/h,"Hz","binding/h","gross spectral limit")
            put(label+"_orbital_period",period,"s","2*pi*hbar/(mass*c^2*alpha^2)","model recurrence, not an independently observed elementary clock")
            put(label+"_gross_3_to_2_vacuum_nm",1/(R*(D(1)/4-D(1)/9))*D("1e9"),"nm","1/[R*(1/4-1/9)]","uncorrected gross vacuum line; not an air wavelength or resolved observed component")
            close(label+": radius ratio",radius/length,1/a)
            close(label+": stationary circular energy derivative",hb*hb/(mass*radius**3),K/radius**2)
            close(label+": total stationary energy",hb*hb/(2*mass*radius**2)-K/radius,-binding)
            close(label+": spectral-energy equivalence",h*c*R,binding)
            close(label+": recurrence/reference-period ratio",period/(h/(mass*c*c)),1/(a*a))
            for n in (2,3,5):
                nn=D(n)
                rn=nn*nn*radius
                en=nn*nn*hb*hb/(2*mass*rn*rn)-K/rn
                close(label+f": n={n} binding progression",en,-binding/(nn*nn))

        # Full equal-end kinetic chart: one diagnostic branch, not fine structure.
        beta = a/2
        root = (1-beta*beta).sqrt()
        momentum = ml*c*beta/root
        r_exact = hb/momentum
        binding_exact = 2*ml*c*c*beta*beta/(1+root)
        binding_quadratic = ml*c*c*a*a/4
        put("equal_end_exact_radius",r_exact,"m","2*hbar*sqrt(1-alpha^2/4)/(m_light*c*alpha)","full two-end kinetic chart, equal masses and n=1, inherited 1/r potential")
        put("equal_end_exact_binding_eV",binding_exact/ev,"eV","2*m_light*c^2*[1-sqrt(1-alpha^2/4)]/eV_to_J","conditional model result, not a fine-structure calculation")
        put("equal_end_binding_fractional_correction",beta*beta/(1+root)**2,"1","binding_exact/binding_quadratic-1","nonlinear kinetic correction within this trial model")
        endpoint_energy=(ml*ml*c**4+c*c*momentum*momentum).sqrt()
        close("full equal-end stationarity",2*c*c*momentum/endpoint_energy,K/hb)
        close("full equal-end stationary total energy",2*endpoint_energy-(K/hb)*momentum,2*ml*c*c*root)
        close("stable exact binding formula",binding_exact,2*ml*c*c*(1-root))
        close("equal-end fractional correction",binding_exact/binding_quadratic-1,rows["equal_end_binding_fractional_correction"]["number"])

        # Back-solving the source norm. This explicitly assumes C_*=K_Q.
        gamma_g = G*ml*mh/(hb*c)
        put("gravitational_pair_dimensionless",gamma_g,"1","G*m_light*m_heavy/(hbar*c)","G and both masses supplied independently of alpha")
        for uu in (-1,0,1):
            u=D(uu)
            rt=(a**4-gamma_g**4*(1-u*u)).sqrt()
            # Stable even when the shift lies below 150 significant digits of a^2.
            delta_x=-gamma_g**4*(1-u*u)/(rt+a*a)-gamma_g*gamma_g*u
            x=rt-gamma_g*gamma_g*u
            frac=delta_x/(a*a)
            base=mu*c*c*a*a/2
            tag={-1:"opposite",0:"perpendicular",1:"aligned"}[uu]
            put("source_Ks_"+tag,mu*c*c*x/2,"J","mu*c^2/2 * [sqrt(alpha^4-gamma_g^4*(1-u^2))-gamma_g^2*u]","conditional reverse circular map C_*=K_Q, n=1; u remains unspecified")
            put("source_Ks_relative_shift_"+tag,frac,"1","(required_Ks - mu*c^2*alpha^2/2)/(mu*c^2*alpha^2/2)","stable small-difference formula; far below physical input precision")
            put("source_Ks_absolute_shift_"+tag,base*frac,"J","base_energy * relative_shift","mathematical correction, not experimentally resolved")
            close("source fourth-power reconstruction, u="+str(uu),gamma_g**4+x*x+2*gamma_g*gamma_g*x*u,a**4)
            close("source reverse through original norm, u="+str(uu),(gamma_g**4/4+x*x/4+gamma_g*gamma_g*x*u/2).sqrt(),a*a/2)
        check("source u=0 shift is negative and retained, not rounded to zero",rows["source_Ks_relative_shift_perpendicular"]["number"]<0)

        # Previous inferred clock; changing zeta is still allowed by these equations.
        Mpl=(hb*c/G).sqrt()
        fpl=Mpl*c*c/h
        tpl=(hb*G/c**5).sqrt()
        kM=G*h*h/c**4
        put("Planck_mass",Mpl,"kg","sqrt(hbar*c/G)","uses measured G; not derived from alpha")
        put("Planck_frequency",fpl,"Hz","M_Pl*c^2/h","uses measured G; cycle frequency, not angular frequency")
        put("light_mass_over_Planck",ml/Mpl,"1","m_light/M_Pl","uses measured mass and G, not a result of alpha")
        put("alpha_G_light",G*ml*ml/(hb*c),"1","G*m_light^2/(hbar*c)","uses mass and G")
        for ztxt in ("0.25","1","4"):
            z=D(ztxt)
            nu=(a/z).sqrt()*fpl
            key=ztxt.replace(".","p")
            put("rate_zeta_"+key,nu,"Hz","sqrt(alpha/zeta)*f_Pl","conditional normalization example; not an independently established completion clock")
            close("rate-product invariance zeta="+ztxt,z*kM*nu*nu,K)
        put("reference_rate_period_in_Planck_times",1/(a.sqrt()*fpl*tpl),"1","2*pi/sqrt(alpha)","requires zeta=1; no count or actual period established")
        close("Planck cycle-frequency convention",fpl*tpl,1/(2*pi))

        # Finite spectral trial has the WRONG SIGN for the declared sense mapping.
        put("spectral_magnitude_product_limit",8*a/pi,"1","8*alpha/pi","magnitude only of (Omega*a/c)*(r/ell); actual required product is negative in this trial")
        for N in (3,10,100):
            nn=D(N)
            angle=pi/(4*nn)
            tan=sine(angle)/cosine(angle)
            req=-2*a/(nn*tan)
            put(f"spectral_required_product_N{N}",req,"1","-2*alpha/[N*tan(pi/(4*N))]","no solution with positive frequency, lengths and c under the original sign assignment")
            close(f"finite spectrum back-substitution N={N}",-nn*req*tan/2,a)
            check(f"finite spectrum requires negative product N={N}",req<0)

        # Printed-input precision and rounding checks, not an independence claim.
        ua=D(INPUTS["alpha"]["standard_uncertainty"])
        put("relative_input_alpha_uncertainty",ua/a,"1","u(alpha)/alpha","measurement precision, irrespective of arithmetic precision")
        put("inverse_alpha_first_order_uncertainty",ua/(a*a),"1","u(alpha)/alpha^2","first-order standard uncertainty from alpha alone")
        put("legacy_N_first_order_uncertainty",4*ua/a**3,"1","4*u(alpha)/alpha^3","first-order standard uncertainty from alpha alone")
        put("inverse_alpha_printed_reference_difference",1/a-D(REFERENCES["inverse_alpha"]),"1","1/(printed alpha) - separately printed inverse alpha","small independent-rounding discrepancy; not physical evidence")
        ref_R=D(REFERENCES["Rydberg_per_m"])
        Rcalc=rows["light_limit_Rydberg"]["number"]
        put("Rydberg_printed_reference_relative_difference",Rcalc/ref_R-1,"1","R_from_printed_mass_and_alpha/R_CODATA-1","correlated adjustment and finite printed rounding; not an independent precision prediction")
        radius_ref=D(REFERENCES["Bohr_radius_m"])
        check("atomic radius agrees with printed NIST value within 1 part in 1e8",abs(rows["light_limit_radius"]["number"]/radius_ref-1)<D("1e-8"))
        check("spectral coefficient agrees within 1 part in 1e8",abs(Rcalc/ref_R-1)<D("1e-8"))
        # Maximum first-order relative change from HALF a last printed digit.
        # This is a rounding allowance, not a propagated standard uncertainty.
        rounding_R=D("5e-42")/ml+2*D("5e-14")/a+D("0.0000005")/ref_R
        put("Rydberg_printed_rounding_allowance",rounding_R,"1","half-last-digit relative allowances for mass + 2*alpha + R","conservative first-order input/output rounding allowance, not statistical uncertainty")
        check("Rydberg small discrepancy is covered by printed decimal rounding",abs(Rcalc/ref_R-1)<rounding_R)
        check("legacy N is separated from its nearest integer beyond alpha uncertainty",abs(2/(a*a)-(2/(a*a)).to_integral_value())>100*rows["legacy_N_first_order_uncertainty"]["number"])
        return rows, checks


def serial(rows):
    return {key: {"value_80_significant_digits": format(val["number"],".79E"),
                  "display_11_significant_digits": format(val["number"],".10E"),
                  **{k:v for k,v in val.items() if k!="number"}}
            for key,val in rows.items()}


def main():
    low, checks = calculate(150)
    high, _ = calculate(230)
    lo, hi = serial(low), serial(high)
    assert lo==hi, [(k,lo[k],hi[k]) for k in lo if lo[k]!=hi[k]]
    checks.append({"name": "all archived 80-significant-digit outputs identical at 150 and 230 working digits", "passed": True})
    provenance = {
        "purpose": "user-requested measured-alpha reverse calculation, not a forward derivation or numerical search",
        "source": SOURCE,
        "source_edition": "2022 CODATA adjustment, NIST table checked 2026-09-15",
        "source_access": "NIST web reader succeeded; a separate urllib download returned HTTP 403, so no full-table snapshot is claimed",
        "input_policy": "the printed central alpha is the sole alpha input; the separately rounded inverse is a reference only",
        "precision_policy": "150 and 230 working Decimal digits, 80 significant output digits; measurement uncertainty retained separately",
        "model_policy": "distinct conditional branches; not imposed simultaneously as one physical theory",
        "python": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "inputs": INPUTS,
        "reference_only_not_inputs": REFERENCES,
    }
    (HERE/"coulomb_calibrated_inputs.json").write_text(json.dumps(provenance,indent=2)+"\n")
    (HERE/"coulomb_calibrated_results.json").write_text(json.dumps({"provenance":provenance,"values":hi,"checks":checks},indent=2)+"\n")
    lines=["# Reverse alpha: 80-significant-digit arithmetic archive", "", "These are computed digits from finite-precision measured inputs, not measured digits. See the [input provenance](./coulomb_calibrated_inputs.json), [reproducible calculation](./coulomb_calibrated.py), and [Coulomb exhibit](../coulomb-from-chirality-directed-depletion.html#research-status) for assumptions and scope.", ""]
    for key,v in hi.items():
        lines += [f"## {key}", "", f"Formula: `{v['formula']}`", "", f"Unit: {v['unit']}. {v['status']}.", "", "```text",v["value_80_significant_digits"],"```",""]
    (HERE/"coulomb_calibrated_precision.md").write_text("\n".join(lines))
    print(f"PASS: {len(checks)} mathematical/input checks; {len(hi)} outputs stable to all 80 archived significant digits.")
    print("No forward alpha derivation or physical model validation is claimed.")
    for key in ("inverse_alpha","sqrt_alpha","network_kappa","K_Q","legacy_N","binding_fraction","connected_phase_degrees","cone_angle_degrees","light_limit_radius","light_limit_binding_eV","light_limit_Rydberg","two_end_radius","two_end_binding_eV","two_end_gross_3_to_2_vacuum_nm","rate_zeta_1","source_Ks_relative_shift_aligned","source_Ks_relative_shift_perpendicular","spectral_magnitude_product_limit","Rydberg_printed_reference_relative_difference","relative_input_alpha_uncertainty"):
        v=hi[key]
        print(key, v["display_11_significant_digits"],v["unit"])


if __name__ == "__main__":
    main()
