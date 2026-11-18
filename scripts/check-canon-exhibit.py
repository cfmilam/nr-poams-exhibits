#!/usr/bin/env python3
"""Fail-closed checks for the generated Canon exhibit and Hub placement."""
from __future__ import annotations

import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import runpy
import subprocess
from urllib.parse import unquote, urlsplit

REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
PAGE = REPO / "the-canon.html"
MANIFEST = REPO / "canon" / "SOURCE-MANIFEST.json"
FORBIDDEN = (
    "Changes from Draft 5", "Effects proposed and not demonstrated",
    "Content owner-ratified on 2026-11-17", "I strongly support the approach",
    "Combined & committed 2026-07-19", "(ratified by Star Lord)",
    "(ratified by Christopher F. Milam)", "Ratified:** 2026-09-24 by Christopher F. Milam",
    "Version receipt", "version receipt", "private owner-review edition",
    "Producer history:", "Verification receipt:", "Quoted corpus-file SHA-256",
    "Integrated substantive Verification passed under", "Owner acceptance, live-context migration and publication remain separate decisions",
    "Draft 3 for Star Lord's approval", "DIRECTION-RESET-2026-11-16.md",
    "owner-ratified public edition", "Owner-approved 15 November 2026",
    "Owner-approved main-session recital", "standard tightened 2026-11-10, owner-ratified",
    "The decisive feature is", "There are three safeguards we must preserve",
    "your ratified constitutive extensions", "all 23 full pages", "final system should generate",
    "The correct final deliverable", "approved and controlled", "approved Project-Start Protocol",
    "The owner's present writing conventions", "this draft has not newly authenticated",
    "later owner prose", "under the owner's amendment", "The owner's thirteenth-dossier amendment",
    "owner prose", "Project-start gate · Draft 3",
)
MAJORS = ["core", "grounding", "canon", "project-start-protocol", "opening-context-block", "canon-notes"]
HEADINGS = ["Recital", "Core", "Grounding", "Canon", "Project Start Protocol", "Opening Context Block", "Canon Notes"]
DOSSIER_KEYS = [f"canon-c{i:02d}" for i in range(1, 14)]
EXPECTED_KEYS = [
    "recital",
    "canon-core",
    "grounding",
    "foundation-wall",
    *DOSSIER_KEYS,
    "project-start-protocol",
    "opening-context",
    "project-context",
    "canon-notes",
]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class Structure(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.h1 = 0
        self.h2s = []
        self.in_h2 = False
        self.links = []
        self.ids: list[str] = []
        self.majors: list[str] = []
        self.cards = 0
        self.refs: list[str] = []
        self.local_anchors: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if tag == "h2":
            self.h2s.append("")
            self.in_h2 = True
        classes = (values.get("class") or "").split()
        if tag == "a":
            self.links.append(values)
        if tag == "h1":
            self.h1 += 1
        if values.get("id"):
            self.ids.append(values["id"] or "")
        if tag == "details" and "major" in (values.get("class") or "").split():
            self.majors.append(values.get("id") or "")
        if tag == "article" and "source-card" in (values.get("class") or "").split():
            self.cards += 1
        for name in ("href", "src"):
            if values.get(name):
                ref = values[name] or ""
                self.refs.append(ref)
                if ref.startswith("#"):
                    self.local_anchors.append(ref[1:])

    def handle_data(self, data):
        if self.in_h2:
            self.h2s[-1] += data

    def handle_endtag(self, tag):
        if tag == "h2":
            self.in_h2 = False


def check_public_texts(page: str, payloads: dict[str, str], manifest: dict) -> None:
    """Check reader payloads and public metadata, not just rendered markup."""
    surfaces = {"page": page, **payloads, "manifest": json.dumps(manifest, ensure_ascii=False)}
    for surface, text in surfaces.items():
        for retired in FORBIDDEN:
            assert retired.casefold() not in text.casefold(), f"{surface}: retired text {retired!r}"


def main() -> None:
    builder = runpy.run_path(str(REPO / "scripts/build-canon-exhibit.py"))
    specs = {spec["key"]: spec for spec in builder["SOURCES"]}
    manifest = json.loads(MANIFEST.read_text())
    assert manifest["schema"] == 2
    assert [row["key"] for row in manifest["sources"]] == EXPECTED_KEYS
    payloads: dict[str, str] = {}
    joined: list[str] = []
    for row in manifest["sources"]:
        snapshot = REPO / "canon" / row["snapshot"]
        data = snapshot.read_bytes()
        assert len(data) == row["bytes"]
        assert sha(data) == row["sha256"]
        payloads[row["key"]] = data.decode()
        joined.append(f"{row['sha256']}  {row['snapshot']}")
        canonical_hash = row.get("canonical_sha256")
        if canonical_hash:
            source = WORKSPACE / row["canonical"]
            source_text = source.read_text()
            assert sha(source.read_bytes()) == canonical_hash
            expected = builder["public_payload"](specs[row["key"]], source_text)
            assert data.decode() == expected
    assert sha(("\n".join(joined) + "\n").encode()) == manifest["source_set_sha256"]
    assert (REPO / "canon/SOURCE-MANIFEST.sha256").read_text() == "\n".join(joined) + "\n"

    page = PAGE.read_text()
    check_public_texts(page, payloads, manifest)
    embedded_match = re.search(
        r'<script id="canon-sources" type="application/json">(.*?)</script>', page, re.S
    )
    assert embedded_match
    assert json.loads(embedded_match.group(1)) == payloads

    structure = Structure()
    structure.feed(page)
    structure.close()
    assert structure.h1 == 1
    assert structure.h2s == HEADINGS, structure.h2s
    assert re.findall(r'<(?:section class="recital"|details class="major") id="([^"]+)"', page) == ["recital"] + MAJORS
    for retired_media in (
        "poams-listen",
        "canon-media-pending",
        "podcast-listen.css",
        "podcast-listen.js",
        'data-podcast="62"',
        "pl-pills",
        "pl-pill",
        "pl-video",
        "pl-player-62",
        "Audio: The Canon",
        "Audio under consideration",
        "Video in NotebookLM: The Canon",
        "https://notebook.google.com/notebook/5174836c-294f-471a-9ebf-f55426c9eb34",
    ):
        assert retired_media not in page, retired_media
    assert not re.search(r"<audio(?:\s|>)", page, re.I)
    for row in manifest["sources"]:
        assert "canon/" + row["snapshot"] in structure.refs
    assert "canon/SOURCE-MANIFEST.json" in structure.refs
    assert sha(payloads["recital"].encode()) == "7c44635f9793612afb6356de06863e972b44d8e580f3a69d8f471471ff00c31e"
    assert payloads["recital"].startswith("Every new exhibit must begin with:")
    assert "activation/migration against this edition is unreviewed" in payloads["canon-notes"]
    recital_html = subprocess.run(["pandoc", str(REPO / "canon/00-recital.md"), "--from=gfm+tex_math_dollars", "--to=html5", "--mathjax", "--shift-heading-level-by=2", "--id-prefix=recital-"], check=True, capture_output=True, text=True).stdout.strip()
    assert '<div class="source-document">' + recital_html + '</div></section>' in page
    assert sum(a.get("href") == "./index.html" and a.get("aria-label") == "Previous: Hub" for a in structure.links) == 2
    assert sum(a.get("href") == "./the-rules-of-persistence.html" and a.get("aria-label") == "Next: The Rules of Persistence" for a in structure.links) == 2
    neighbor = Structure()
    neighbor.feed((REPO / "the-rules-of-persistence.html").read_text())
    assert any(a.get("href") == "./the-canon.html" and a.get("aria-label") == "Previous: The Canon" for a in neighbor.links)
    assert "./the-world-is-language.html" not in neighbor.refs
    assert "./the-designers-reading.html" in neighbor.refs
    assert structure.majors == MAJORS
    assert structure.cards == 20
    assert len(structure.ids) == len(set(structure.ids))
    assert {unquote(a) for a in structure.local_anchors} <= set(structure.ids)
    for ref in structure.refs:
        parsed = urlsplit(ref)
        if parsed.scheme or parsed.netloc or ref.startswith("#"):
            continue
        target = (REPO / unquote(parsed.path)).resolve()
        assert target.exists(), f"missing local dependency: {ref}"
        if parsed.fragment and target.suffix == ".html":
            linked = Structure()
            linked.feed(target.read_text())
            assert unquote(parsed.fragment) in linked.ids, ref

    scripts = []
    for attrs, body in re.findall(r"<script([^>]*)>(.*?)</script>", page, re.S):
        if "src=" not in attrs and "application/json" not in attrs:
            scripts.append(body)
    subprocess.run(["node", "--check", "-"], input="\n".join(scripts), text=True, check=True)

    hub = (REPO / "index.html").read_text()
    foundation_positions = [
        hub.index(">The Logos<"),
        hub.index(">The World Is Language<"),
        hub.index(">The Rules of Persistence<"),
        hub.index(">POAMS Historical Framework<"),
        hub.index(">The Canon<"),
    ]
    assert foundation_positions == sorted(foundation_positions)
    hub_structure = Structure()
    hub_structure.feed(hub)
    pages = {
        ref.split("#", 1)[0].split("?", 1)[0]
        for ref in hub_structure.refs
        if ref.startswith("./")
        and ref.split("#", 1)[0].split("?", 1)[0].endswith(".html")
        and ref.split("#", 1)[0].split("?", 1)[0] != "./index.html"
    }
    fallback = int(re.search(r'id="exhibit-count">(\d+)<', hub).group(1))
    assert fallback == len(pages)
    assert "./the-canon.html" in pages

    before = PAGE.read_bytes()
    subprocess.run(["python3", str(REPO / "scripts" / "build-canon-exhibit.py")], check=True)
    assert PAGE.read_bytes() == before
    print(
        json.dumps(
            {
                "status": "pass",
                "source_set_sha256": manifest["source_set_sha256"],
                "sources": len(payloads),
                "source_cards": structure.cards,
                "unique_ids": len(structure.ids),
                "hub_linked_pages": len(pages),
                "deterministic_build": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
