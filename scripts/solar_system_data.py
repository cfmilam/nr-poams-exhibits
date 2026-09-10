#!/usr/bin/env python3
"""Build the JPL measurement snapshot used by ``solar-system.html``.

This script does no resonance fitting.  It fixes one observable for every entry:
the sidereal orbital period (the period of orbital writhe).  The resulting JSON
also carries the apsidal and nodal periods published in JPL's satellite mean-
elements table so that a physical resonant argument can include every angular
rate it actually uses.

Run from anywhere::

    python3 poams-exhibits-public/scripts/solar_system_data.py

The output is deterministic apart from JPL source updates and the retrieval date.
"""

from __future__ import annotations

import html
import json
import re
from datetime import date
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "solar-system-jpl.json"
JS_OUTPUT = ROOT / "solar-system-jpl.js"

SOURCES = {
    "planet_physical": "https://ssd.jpl.nasa.gov/planets/phys_par.html",
    "satellite_physical": "https://ssd.jpl.nasa.gov/sats/phys_par/",
    "satellite_elements": "https://ssd.jpl.nasa.gov/sats/elem/sep.html",
}


class TableParser(HTMLParser):
    """Small tolerant table parser for JPL's intentionally non-strict HTML."""

    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self.table: list[list[str]] | None = None
        self.row: list[str] | None = None
        self.cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        del attrs
        if tag == "table" and self.table is None:
            self.table = []
        elif self.table is not None and tag == "tr":
            self.row = []
        elif self.table is not None and tag in ("td", "th"):
            # Several JPL cells omit </td>; a new cell therefore closes the old one.
            self._finish_cell()
            self.cell = []
        elif self.cell is not None and tag == "br":
            self.cell.append(" ")

    def handle_endtag(self, tag: str) -> None:
        if self.table is not None and tag in ("td", "th"):
            self._finish_cell()
        elif self.table is not None and tag == "tr":
            self._finish_cell()
            if self.row:
                self.table.append(self.row)
            self.row = None
        elif tag == "table" and self.table is not None:
            self.tables.append(self.table)
            self.table = None

    def handle_data(self, data: str) -> None:
        if self.cell is not None:
            self.cell.append(data)

    def _finish_cell(self) -> None:
        if self.cell is not None and self.row is not None:
            self.row.append(" ".join(html.unescape("".join(self.cell)).split()))
        self.cell = None


def fetch(url: str) -> str:
    request = Request(url, headers={"User-Agent": "POAMS exhibit data audit/1.0"})
    with urlopen(request, timeout=45) as response:
        return response.read().decode("utf-8")


def tables(document: str) -> list[list[list[str]]]:
    parser = TableParser()
    parser.feed(document)
    return parser.tables


def number(value: str) -> float | None:
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", value.replace(",", ""))
    return float(match.group()) if match else None


def planet_rows(document: str) -> list[dict[str, object]]:
    parsed = tables(document)
    rows: list[dict[str, object]] = []
    for family, table in zip(("planet", "dwarf planet"), parsed[:2], strict=True):
        for row in table[2:]:
            if len(row) < 7:
                continue
            rows.append(
                {
                    "host": "Sun",
                    "name": row[0],
                    "family": family,
                    "radius_km": number(row[2]),
                    "mass_table_units": number(row[3]),
                    "rotation_days": number(row[5]),
                    "period_days": number(row[6]) * 365.25,
                }
            )
    return rows


def satellite_rows(physical_document: str, elements_document: str) -> list[dict[str, object]]:
    physical = tables(physical_document)[0]
    by_code: dict[str, dict[str, object]] = {}
    for row in physical[2:]:
        if len(row) < 10:
            continue
        by_code[row[2]] = {
            "host": row[0],
            "name": row[1],
            "code": row[2],
            "family": "satellite",
            "gm_km3_s2": number(row[3]),
            "radius_km": number(row[6]),
        }

    joined: list[dict[str, object]] = []
    seen: set[str] = set()
    # The final table is the legend.  A code can occur in only one active solution
    # among the characterized satellites; first occurrence is the main table.
    for table in tables(elements_document)[:-1]:
        for row in table[1:]:
            if len(row) < 11 or row[1] not in by_code or row[1] in seen:
                continue
            period = number(row[8])
            if period is None:
                continue
            seen.add(row[1])
            joined.append(
                {
                    **by_code[row[1]],
                    "a_km": number(row[2]),
                    "eccentricity": number(row[3]),
                    "inclination_deg": number(row[6]),
                    "period_days": abs(period),
                    "apsis_period_years": number(row[9]),
                    "node_period_years": number(row[10]),
                    # UI-only operational cut; every characterized satellite remains
                    # in the snapshot and can be shown regardless of this label.
                    "prominence": "major" if (by_code[row[1]]["radius_km"] or 0) >= 190 else "characterized",
                }
            )

    missing = sorted(set(by_code) - seen)
    if missing:
        raise RuntimeError(f"satellite element rows missing for codes: {missing}")
    return joined


def main() -> None:
    documents = {name: fetch(url) for name, url in SOURCES.items()}
    planets = planet_rows(documents["planet_physical"])
    satellites = satellite_rows(documents["satellite_physical"], documents["satellite_elements"])

    assert len(planets) == 13, len(planets)
    assert len([body for body in planets if body["family"] == "planet"]) == 8
    assert len(satellites) >= 46, len(satellites)
    assert {"Moon", "Io", "Titan", "Triton", "Charon"} <= {body["name"] for body in satellites}

    payload = {
        "schema": "poams.solar-system-jpl.v1",
        "retrieved": date.today().isoformat(),
        "observable": "sidereal orbital period / orbital angular rate",
        "source_notes": [
            "JPL mean satellite elements describe fitted precessing ellipses and are not ephemerides.",
            "Near an integer period ratio is a commensurability screen, not proof of resonant-angle libration.",
            "A physical resonance verdict must include the relevant apsidal/nodal rates and a channel-width or libration test.",
        ],
        "sources": SOURCES,
        "heliocentric": sorted(planets, key=lambda body: body["period_days"]),
        "satellites": sorted(satellites, key=lambda body: (body["host"], body["period_days"])),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    JS_OUTPUT.write_text(
        "// Generated by scripts/solar_system_data.py; do not edit by hand.\n"
        + "window.SOLAR_SYSTEM_JPL = "
        + json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )

    major_count = sum(body["prominence"] == "major" for body in satellites)
    print(f"wrote {OUTPUT}")
    print(f"wrote {JS_OUTPUT}")
    print(f"heliocentric bodies: {len(planets)} (8 planets, 5 dwarf planets)")
    print(f"characterized satellites: {len(satellites)} ({major_count} radius >= 190 km)")


if __name__ == "__main__":
    main()
