#!/usr/bin/env python3
"""Independent arithmetic and source guards for the Rate of Time exhibit."""

from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 50

C = Decimal("299792458")
DAY = Decimal("86400")
GM_E = Decimal("3.986004418e14")
R_E = Decimal("6378137")
V_SPIN = Decimal("465.1013")


def near(label: str, actual: Decimal, expected: Decimal, tolerance: Decimal) -> None:
    error = abs(actual - expected)
    if error > tolerance:
        raise AssertionError(f"{label}: {actual} != {expected} (error {error})")
    print(f"PASS {label}: {actual}")


beta = Decimal("0.99942")
rate = (Decimal(1) - beta * beta).sqrt()
gamma = Decimal(1) / rate
near("muon gamma", gamma, Decimal("29.37"), Decimal("0.01"))
near("muon laboratory mean life / us", Decimal("2.1970") * gamma,
     Decimal("64.5"), Decimal("0.1"))
near("muon turns per mean life", Decimal("2.1970") * gamma / Decimal("0.149"),
     Decimal("433"), Decimal("1"))

rotor = Decimal(340) ** 2 / (Decimal(2) * C * C)
near("rotor deficit", rotor, Decimal("6.43e-13"), Decimal("0.01e-13"))

jet_motion = Decimal(250) ** 2 / (Decimal(2) * C * C)
jet_altitude = Decimal("9.80665") * Decimal(11000) / (C * C)
jet_ns_day = (jet_altitude - jet_motion) * DAY * Decimal("1e9")
near("jet example / ns per day", jet_ns_day, Decimal("74"), Decimal("1"))

gps_r = Decimal("26562000")
gps_orbit = Decimal(3) * GM_E / (Decimal(2) * gps_r * C * C)
ground = GM_E / (R_E * C * C) + V_SPIN ** 2 / (Decimal(2) * C * C)
gps_delta = ground - gps_orbit
near("GPS orbit deficit", gps_orbit, Decimal("2.50e-10"), Decimal("0.01e-10"))
near("ground deficit", ground, Decimal("6.97e-10"), Decimal("0.01e-10"))
near("GPS net / us per day", gps_delta * DAY * Decimal("1e6"),
     Decimal("38.6"), Decimal("0.1"))

iss_r = R_E + Decimal(420000)
iss_delta = ground - Decimal(3) * GM_E / (Decimal(2) * iss_r * C * C)
near("ISS example / us per day", iss_delta * DAY * Decimal("1e6"),
     Decimal("-24"), Decimal("1"))

lhc_rate = Decimal(1) / Decimal(6930)
span_ratio = (Decimal(1) - lhc_rate) / rotor
near("common-deficit span ratio", span_ratio, Decimal("1.55e12"), Decimal("0.01e12"))

factory_mhz = Decimal("10.23") * (Decimal(1) - Decimal("4.47e-10"))
near("GPS rounded factory frequency / MHz", factory_mhz,
     Decimal("10.22999999543"), Decimal("0.00000000001"))

source = Path(__file__).resolve().parents[1] / "the-rate-of-time.html"
html = source.read_text(encoding="utf-8")
required = [
    'max="99.99" step="0.001"',
    "bs.value='99.942'",
    "const C=299792458, maxBeta=0.9999",
    "Inputs, computed outputs, specifications",
    "about twelve orders of magnitude",
    "three halves as large as",
    "Prophetic horizon:",
]
for token in required:
    if token not in html:
        raise AssertionError(f"missing source guard: {token}")
if "Math.min(v/3e8" in html or "Math.pow(10,p*9)" in html:
    raise AssertionError("forbidden superluminal-cap implementation remains")
if "every row measured" in html or "factor of three</strong> over the static-hover" in html:
    raise AssertionError("known evidence or comparison overclaim remains")

print("PASS source guards")
print("RATE OF TIME AUDIT PASS")
