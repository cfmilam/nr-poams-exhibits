#!/usr/bin/env python3
"""Deterministic arithmetic/source audit for The Winding Register."""

from __future__ import annotations

import cmath
import math
from html.parser import HTMLParser
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "the-winding-register.html"


class PageParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: list[str] = []
        self.local_links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if values.get("id"):
            self.ids.append(values["id"] or "")
        href = values.get("href")
        if tag == "a" and href and not href.startswith(("http://", "https://", "#", "mailto:")):
            self.local_links.append(href)


def q_tail(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def main() -> None:
    source = PAGE.read_text(encoding="utf-8")

    # A modulo phase cannot identify the completed-turn count.
    phi = 0.731
    for n in range(-20, 21):
        assert abs(cmath.exp(1j * (phi + 2 * math.pi * n)) - cmath.exp(1j * phi)) < 2e-14

    # Reproduce the page's raw nearest-integer Gaussian error formula.
    sigma = 0.1
    p_raw = 2.0 * q_tail(0.5 / sigma)
    assert math.isclose(p_raw, 5.73303143758389e-7, rel_tol=2e-14)
    commits = 100_000
    p_fail = 1.0 - (1.0 - p_raw) ** commits
    assert 0.055 < p_fail < 0.056
    assert abs(p_fail - commits * p_raw) < 0.002

    # Carry is a radix-boundary event; it does not occur on every unit increment.
    radix = 10
    carries = []
    for value in range(30):
        low = value % radix
        next_low = (value + 1) % radix
        if low == radix - 1 and next_low == 0:
            carries.append(value + 1)
    assert carries == [10, 20, 30]

    required = (
        "Integer winding gives $\\mathbb Z$, not automatically $\\mathbb Z/N\\mathbb Z$",
        "Interference returns a residue, not an unbounded sum",
        "Under independent raw commit errors",
        "correlated drift requires its own model",
        "proposed architecture and finite test program",
        "doi:10.1109/JRPROC.1959.287195",
        "doi:10.1088/1367-2630/9/10/370",
    )
    for text in required:
        assert text in source, text

    parser = PageParser()
    parser.feed(source)
    assert len(parser.ids) == len(set(parser.ids)), "duplicate HTML id"
    for href in parser.local_links:
        target = (ROOT / href.split("#", 1)[0]).resolve()
