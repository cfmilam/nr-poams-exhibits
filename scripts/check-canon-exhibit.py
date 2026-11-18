#!/usr/bin/env python3
"""Fail-closed checks for the generated Canon exhibit and Hub placement."""
from __future__ import annotations

import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit

REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
PAGE = REPO / "the-canon.html"
MANIFEST = REPO / "canon" / "SOURCE-MANIFEST.json"
FORBIDDEN = ("Changes from Draft 5", "Effects proposed and not demonstrated")
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
        self.pending = 0
        self.players = 0
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
        self.pending += int("canon-media-pending" in classes)
        self.players += int("poams-listen" in classes)
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


def public_protocol(text: str) -> str:
    marker = "\n## The fiber program, replayed\n"
    assert marker in text
    return text.split(marker, 1)[0].rstrip() + "\n"


def public_opening(text: str) -> str:
    return "".join(
        line
        for line in text.splitlines(keepends=True)
        if not line.startswith("<!-- GROUNDING:")
        and not line.startswith("**Current work and its state:**")
    ).strip() + "\n"


def main() -> None:
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
            expected = source_text
            if row["key"] == "project-start-protocol":
                expected = public_protocol(source_text)
            elif row["key"] == "opening-context":
                expected = public_opening(source_text)
            assert data.decode() == expected
    assert sha(("\n".join(joined) + "\n").encode()) == manifest["source_set_sha256"]
    assert (REPO / "canon/SOURCE-MANIFEST.sha256").read_text() == "\n".join(joined) + "\n"

    page = PAGE.read_text()
    for heading in FORBIDDEN:
        assert heading not in page
        assert all(heading not in text for text in payloads.values())
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
    if not (REPO / "canon/audio-binding.json").exists():
        assert structure.pending == 1 and structure.players == 0
        assert not re.search(r'<script[^>]+src="[^"]*podcast-listen.js', page)
        assert page.count('aria-disabled="true" title=') == 2
    else:
        assert structure.players == 1 and structure.pending == 0
        assert '<script defer src="./podcast-listen.js"></script>' in page
    for row in manifest["sources"]:
        assert "canon/" + row["snapshot"] in structure.refs
    assert "canon/SOURCE-MANIFEST.json" in structure.refs
    assert sha(payloads["recital"].encode()) == "7054199e84e95df252e4e6eb7f32e41bd10d370288be7ad364c35d4e25c03df7"
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
    assert 'data-podcast="62"' in page
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
    assert hub.index(">The Canon<") < hub.index(">The Rules of Persistence<")
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
