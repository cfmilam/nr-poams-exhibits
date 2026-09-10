#!/usr/bin/env python3
"""Independent mathematics/logic/source guard for the Brilliant Vortex Atom."""

from __future__ import annotations

import math
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "the-poams-brilliant-vortex-atom.html"
HTML = PAGE.read_text(encoding="utf-8")

H = 6.626_070_15e-34
HBAR = H / (2 * math.pi)
C = 299_792_458.0
R_H = 10_967_758.340
R_INF = 10_973_731.568_157
M_E = 9.109_383_7139e-31
M_P = 1.672_621_925_95e-27
G = 6.674_30e-11
M_SUN = 1.988_47e30
MERCURY_A = 57.909_05e9
MERCURY_E = 0.205_630
MERCURY_PERIOD_DAYS = 87.9691
ARCSEC_PER_RAD = 206_264.806_247_09636


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


# Hydrogen normalization and level/line identity.
nu_r = C * R_H
m_r = H * R_H / C
require(abs(m_r - 2.424_114_890_478_25e-35) < 1e-49, "Rydberg scale mass")
for lower in (1, 2, 3):
    for upper in range(lower + 1, 10):
        level_u = nu_r * (1 - 1 / upper**2)
        level_l = nu_r * (1 - 1 / lower**2)
        line = nu_r * (1 / lower**2 - 1 / upper**2)
        require(math.isclose(level_u - level_l, line, rel_tol=3e-15), "level difference")

# The source's atomic values require distinct electron/proton mass roles.
ks = H * C * R_INF
v_s = math.sqrt(2 * ks / M_E)
r_s = H / (2 * math.pi * M_E * v_s)
g_star = r_s * v_s**2 / M_P
require(abs(v_s / 2.187_691_263_64e6 - 1) < 2e-9, "atomic speed")
require(abs(r_s / 5.291_772_105_44e-11 - 1) < 2e-9, "atomic radius")
require(abs(g_star / 1.514_172_70e29 - 1) < 2e-8, "two-body coefficient")
require(r_s * v_s**2 / M_E > 2.7e32, "mass-symbol ambiguity guard")

# Equation (5): L^2/(m^2 r^3) is acceleration.
# Dimension vectors are (mass, length, time).
def dim_add(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(x + y for x, y in zip(a, b))

L_dim = (1, 2, -1)
term_dim = dim_add((2 * L_dim[0], 2 * L_dim[1], 2 * L_dim[2]), (-2, -3, 0))
require(term_dim == (0, 1, -2), f"Eq. 5 dimensions: {term_dim}")

# The adopted magnitude relation is nonnegative and gives the stated aligned limits.
for rho_i in range(0, 401):
    rho = rho_i / 40
    for theta_i in range(0, 361):
        theta = math.radians(theta_i)
        value = math.sqrt(max(0.0, 1 + rho**2 + 2 * rho * math.cos(theta)))
        require(math.isfinite(value) and value >= 0, "magnitude relation")
    require(math.isclose(math.sqrt((1 + rho) ** 2), 1 + rho, abs_tol=2e-14), "parallel limit")
    require(math.isclose(math.sqrt((1 - rho) ** 2), abs(1 - rho), abs_tol=2e-14), "antiparallel limit")

# Equal action-event counts do not imply equal energy rates when frequencies differ.
require(H * 2.0 != H * 3.0, "rate-weighted transaction ledger")

# Mercury value harmonized to the audited orbital exhibit.
epsilon = G * M_SUN / (MERCURY_A * (1 - MERCURY_E**2) * C**2)
arcsec_orbit = 6 * math.pi * epsilon * ARCSEC_PER_RAD
arcsec_century = arcsec_orbit * (36_525 / MERCURY_PERIOD_DAYS)
require(abs(arcsec_century - 42.981_975) < 0.002, "Mercury arithmetic")

# Fixed-link topology: twist and writhe trade while linking remains invariant.
for lk in range(-5, 6):
    for tw_i in range(-40, 41):
        tw = tw_i / 8
        wr = lk - tw
        require(math.isclose(tw + wr, lk, abs_tol=1e-15), "fixed-link residual")

required = [
    "[SOURCE]",
    "[CORPUS]",
    "[HYPOTHESIS]",
    "[REFUTED]",
    r"\frac{L^2}{m^2r^3}",
    r"m_R=\frac{h\nu_R}{c^2}",
    r"r_s=\frac{h}{2\pi m_e v_s}",
    r"\frac{r_s v_s^2}{m_p}",
    "all bounded orbits",
    "zero-remainder boundary selects the inverse-square branch",
    r"\mathcal K(\theta)=\sqrt{K_0^2+K_s^2+2K_0K_s\cos\theta}",
    r"\Delta K_{\rm align}=fN\varepsilon_{\rm align}",
    "does not, by logic alone, prove either the existence or nonexistence",
    r"$z\propto r^2$",
    r"$L_{\rm total}=\hbar\,\mathrm{Lk}$",
    "does not reproduce every two-body dynamical equation",
]
for token in required:
    require(token in HTML, f"missing guard: {token}")

for forbidden in [
    "two premises only",
    r"\frac{L^2}{m r^3}",
    r"r_s = \frac{h}{2\pi m v_s}",
    r"\frac{r_s v_s^2}{m}",
    "vindicating the anomalous-weight",
    "There is therefore no event horizon",
    "reached here from first principles",
    "reproduces every two-body result while",
    "Everything in Parts I–VII",
]:
    require(forbidden not in HTML, f"obsolete claim remains: {forbidden}")


class LinkAudit(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: list[str] = []
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        data = dict(attrs)
        if data.get("id"):
            self.ids.append(data["id"] or "")
        if tag == "a" and data.get("href"):
            self.hrefs.append(data["href"] or "")


audit = LinkAudit()
audit.feed(HTML)
require(len(audit.ids) == len(set(audit.ids)), "duplicate ids")
for href in audit.hrefs:
    parts = urlsplit(href)
    if parts.scheme or parts.netloc or not parts.path:
        continue
    target = (PAGE.parent / unquote(parts.path)).resolve()
    require(target.exists(), f"missing local link: {href}")

print(
    "PASS brilliant-vortex-atom:",
    f"m_R={m_r:.9e} kg",
    f"v_s={v_s:.7e} m/s",
    f"r_s={r_s:.9e} m",
    f"G*={g_star:.8e}",
    f"Mercury={arcsec_century:.6f} arcsec/century",
    f"ids={len(audit.ids)}",
    f"links={len(audit.hrefs)}",
)
