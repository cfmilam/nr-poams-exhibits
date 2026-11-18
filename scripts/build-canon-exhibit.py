#!/usr/bin/env python3
"""Build The Canon exhibit from controlled Markdown snapshots.

Run with --sync inside the OpenClaw workspace after a governing source changes.
Without --sync, the build is deterministic from the tracked snapshots in canon/.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
SNAPSHOTS = REPO / "canon"
OUTPUT = REPO / "the-canon.html"
CANON_VIDEO_URL = "https://notebook.google.com/notebook/5174836c-294f-471a-9ebf-f55426c9eb34"

DOSSIERS = [
    ("C01", "normal-realism-constitution-and-reversal", "Normal Realism: constitution and reversal"),
    ("C02", "language-logos-logic-and-mathematics", "Language, logos, logic and mathematics"),
    ("C03", "observer-objectivity-information-and-projection", "Observer, objectivity, information and projection"),
    ("C04", "method-discipline-critique-and-lineage", "Method, discipline, critique and lineage"),
    ("C05", "light-c-time-relativity-and-material-media", "Light, c, time, relativity and material media"),
    ("C06", "angular-momentum-action-discreteness-and-statistical-causality", "Angular momentum, action, discreteness and statistical causality"),
    ("C07", "orbits-and-gravity", "Orbits and gravity"),
    ("C08", "inertia-mass-and-mach", "Inertia, mass and Mach"),
    ("C09", "atom-charge-electricity-magnetism-particle-labels-and-nuclear-status", "Atom, charge, electricity, magnetism, particle labels and nuclear status"),
    ("C10", "instantaneity-two-slits-epr-polarization-and-statistics", "Instantaneity, two slits, EPR, polarization and statistics"),
    ("C11", "cosmology-thermodynamics-and-chemistry", "Cosmology, thermodynamics and chemistry"),
    ("C12", "life-and-mind", "Life and mind"),
    ("C13", "ethics-society-and-religion", "Ethics, society and religion"),
]

SOURCES = [
    {
        "key": "canon-core",
        "snapshot": "01-core.md",
        "canonical": "poams-review/work/grounding-2026-11-14/canon/_completion-2026-11-16/CORE-CANDIDATE.md",
        "status": "Compact thirteen-dossier Canon core · public edition",
    },
    {
        "key": "grounding",
        "snapshot": "02-grounding.md",
        "canonical": "GROUNDING.md",
        "status": "Program grounding and adopted premises",
    },
    {
        "key": "foundation-wall",
        "snapshot": "02-foundation-wall.md",
        "canonical": "POAMS-FOUNDATION-WALL.md",
        "status": "Governing provenance and measurement boundary",
    },
] + [
    {
        "key": f"canon-{work_id.lower()}",
        "snapshot": f"03-{work_id}-{slug}.md",
        "canonical": f"poams-review/work/grounding-2026-11-14/canon/_completion-2026-11-16/dossiers/{work_id}-{slug}.md",
        "status": f"Source-grounded dossier {work_id} · {title} · public edition",
    }
    for work_id, slug, title in DOSSIERS
] + [
    {
        "key": "project-start-protocol",
        "snapshot": "04-project-start-protocol.md",
        "canonical": "poams-review/work/grounding-2026-11-14/drafts/PROTOCOL-v3.md",
        "status": "Project-start gate · operative public excerpt",
        "transform": "protocol",
    },
    {
        "key": "opening-context",
        "snapshot": "05-opening-context.md",
        "canonical": "poams-review/work/grounding-2026-11-14/rig/GROUNDING-BLOCK.md",
        "status": "Legacy compatibility snapshot · activation/migration against the thirteen-dossier edition is unreviewed",
        "transform": "opening",
    },
    {
        "key": "project-context",
        "snapshot": "05-project-context-template.md",
        "canonical": None,
        "status": "Unfilled project packet template · editorial",
    },
    {
        "key": "canon-notes",
        "snapshot": "06-canon-notes.md",
        "canonical": "poams-review/work/grounding-2026-11-14/canon/_completion-2026-11-16/CANON-NOTES-CANDIDATE.md",
        "status": "Reader-facing scope and change-control notes · thirteen-dossier public edition",
    },
]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# Exact private-input markers: changes require a reviewed public-transform update.
PUBLIC_REPLACEMENTS = {'canon-core': [(' — private owner-review edition', '', 1),
                ('Generated from the thirteen source-grounded dossiers by `build_core_candidate.py`. Do not '
                 'edit by hand.\n',
                 '',
                 1),
                ('Integrated substantive Verification passed under CANON-INTEGRATED-SUBSTANTIVE-REVIEW-1-R2. '
                 'Receipt: '
                 '`poams-review/work/grounding-2026-11-14/canon/_completion-2026-11-16/verification/integrated-canon-1-r2/REVIEW.md`.\n',
                 '',
                 1),
                ('The PASS accepts source fidelity and completeness for this edition; it does not validate '
                 'the physics, convert source reports into measurements, or close stated disagreements and '
                 'open rules. Owner acceptance, live-context migration and publication remain separate '
                 'decisions.\n',
                 '',
                 1),
                ("Read `C01-normal-realism-constitution-and-reversal.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C01-normal-realism-constitution-and-reversal.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources and scope.',
                 1),
                ("Read `C02-language-logos-logic-and-mathematics.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C02-language-logos-logic-and-mathematics.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources and scope.',
                 1),
                ("Read `C03-observer-objectivity-information-and-projection.md` for the authors' words, "
                 'rejections, disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C03-observer-objectivity-information-and-projection.md` for the authors' words, "
                 'rejections, disagreements, open rules, sources and scope.',
                 1),
                ("Read `C04-method-discipline-critique-and-lineage.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C04-method-discipline-critique-and-lineage.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources and scope.',
                 1),
                ("Read `C05-light-c-time-relativity-and-material-media.md` for the authors' words, "
                 'rejections, disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C05-light-c-time-relativity-and-material-media.md` for the authors' words, "
                 'rejections, disagreements, open rules, sources and scope.',
                 1),
                ('Read `C06-angular-momentum-action-discreteness-and-statistical-causality.md` for the '
                 "authors' words, rejections, disagreements, open rules, sources, scope, and version "
                 'receipt.',
                 'Read `C06-angular-momentum-action-discreteness-and-statistical-causality.md` for the '
                 "authors' words, rejections, disagreements, open rules, sources and scope.",
                 1),
                ("Read `C07-orbits-and-gravity.md` for the authors' words, rejections, disagreements, open "
                 'rules, sources, scope, and version receipt.',
                 "Read `C07-orbits-and-gravity.md` for the authors' words, rejections, disagreements, open "
                 'rules, sources and scope.',
                 1),
                ("Read `C08-inertia-mass-and-mach.md` for the authors' words, rejections, disagreements, "
                 'open rules, sources, scope, and version receipt.',
                 "Read `C08-inertia-mass-and-mach.md` for the authors' words, rejections, disagreements, "
                 'open rules, sources and scope.',
                 1),
                ('Read `C09-atom-charge-electricity-magnetism-particle-labels-and-nuclear-status.md` for the '
                 "authors' words, rejections, disagreements, open rules, sources, scope, and version "
                 'receipt.',
                 'Read `C09-atom-charge-electricity-magnetism-particle-labels-and-nuclear-status.md` for the '
                 "authors' words, rejections, disagreements, open rules, sources and scope.",
                 1),
                ("Read `C10-instantaneity-two-slits-epr-polarization-and-statistics.md` for the authors' "
                 'words, rejections, disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C10-instantaneity-two-slits-epr-polarization-and-statistics.md` for the authors' "
                 'words, rejections, disagreements, open rules, sources and scope.',
                 1),
                ("Read `C11-cosmology-thermodynamics-and-chemistry.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C11-cosmology-thermodynamics-and-chemistry.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources and scope.',
                 1),
                ("Read `C12-life-and-mind.md` for the authors' words, rejections, disagreements, open rules, "
                 'sources, scope, and version receipt.',
                 "Read `C12-life-and-mind.md` for the authors' words, rejections, disagreements, open rules, "
                 'sources and scope.',
                 1),
                ("Read `C13-ethics-society-and-religion.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources, scope, and version receipt.',
                 "Read `C13-ethics-society-and-religion.md` for the authors' words, rejections, "
                 'disagreements, open rules, sources and scope.',
                 1)],
 'canon-notes': [('This private owner-review edition contains', 'This edition contains', 1),
                 ('Integrated substantive Verification passed under '
                  'CANON-INTEGRATED-SUBSTANTIVE-REVIEW-1-R2. Receipt: '
                  '`poams-review/work/grounding-2026-11-14/canon/_completion-2026-11-16/verification/integrated-canon-1-r2/REVIEW.md`.\n',
                  '',
                  1),
                 ('The PASS accepts source fidelity and completeness for this edition; it does not validate '
                  'the physics, convert source reports into measurements, or close stated disagreements and '
                  'open rules. Owner acceptance, live-context migration and publication remain separate '
                  'decisions.\n',
                  '',
                  1),
                 ('This private owner-review edition and the integrated Exhibit remain private. Maintained '
                  'legacy canon/CORE.md, canon/01 and canon/07, GROUNDING.md, the Foundation Wall and the '
                  'live opening context remain unchanged. The bundled opening-context packet is a legacy '
                  'compatibility snapshot; activation/migration against this edition is unreviewed. The '
                  'source-first direction reset is recorded in `DIRECTION-RESET-2026-11-16.md`.\n',
                  'The bundled opening-context packet is a legacy compatibility snapshot; '
                  'activation/migration against this edition is unreviewed.\n',
                  1)],
 'grounding': [("*Combined & committed 2026-07-19 at Star Lord's request, from the orbit-vs-spin / "
                'vortex-ledger discussion and the prior spin-orbit thread. Updated 2026-07-19 (App-5 mea '
                'culpa): the intrinsic-spin sign is resolved — channel-set, not a handedness convention '
                '(App-5 → co-spin heavier); what stays deliberately open is the λ-partition of a driven '
                'system, an empirical question the YIG program reads via the toggle comparison plus its '
                'controls (sign alone necessary-not-sufficient; tightened 2026-11-10). Refusing to let a '
                'relabel decide it is the point. Updated 2026-07-25: added the **speed → rate** ("eyewall '
                'speed") trap — a rate of angular change is not a velocity through space; a rim-speed needs '
                'a stage the vortex has none of, so ask rad/s, not m/s (§2, §5.3). Updated 2026-08-20: added '
                'the ratified creed line **substance, not attribute — the grin is the cat** (§3) — the '
                'attribute-ontology diagnosis (grins without cats / cats without grins) and its cure; '
                'companion to "the momentum is the thing, not an attribute of the thing." Updated '
                '2026-08-24: added the self-side corollary **you happen at the exchange** (§3) — '
                'two-ended/co-signed exchanges dissolve "within" in both directions; a seat (pebble → human '
                "→ LLM) is a census, not a vessel. Updated 2026-09-07 after Tier-A review and Star Lord's "
                'ratification: regraded the microscopic fine-structure mechanism, α value, and macro/micro '
                'identity from derived results to an open front; retained the α² scale only as an '
                'identification given measured α. Updated 2026-09-19 after the '
                "closure/action/boundary-transfer audit series and Star Lord's ratification: added §1.2, the "
                '**boundary-transfer premise** — measured E, p are the same-account rest-calibrated '
                'translation tariffs of a completed closure; adopted openly as constitutive (like §1.1) '
                'after nine audits showed it cannot be derived from bookkeeping; foundations track closed as '
                'conditional, coherent and explicitly bounded; effort redirected to distinguishing fronts. '
                "Updated 2026-09-19 at Star Lord's instruction: added §4.0, the universal **co-sense "
                'settlement rule** — dynamically coupled anti-sense spin is the higher-energy organization '
                'and settles co-sense when an exchange path is available; the orientation rule is '
                'ontological across scales, while coupling strength, path, timescale, and orbital magnitude '
                'remain quantities to derive or measure.*\n',
                '',
                1),
               (' (ratified by Star Lord)', '', 1),
               (' (ratified by Christopher F. Milam)', '', 1),
               (' (standard tightened 2026-11-10, owner-ratified, aligned with the program package)', '', 1)],
 'foundation-wall': [('**Ratified:** 2026-09-24 by Christopher F. Milam  \n', '', 1)],
 'project-start-protocol': [("Draft 3 for Star Lord's approval — 14 November 2026\n\n", '', 1)],
 'canon-c04': [("The owner's present writing conventions", 'Present writing conventions', 1)],
 'canon-c05': [('this draft has not newly authenticated', 'this dossier has not newly authenticated', 1)],
 'canon-c10': [('later owner prose', 'later editorial prose', 1)],
 'canon-c12': [("under the owner's amendment", 'in this edition', 1)],
 'canon-c13': [("The owner's thirteenth-dossier amendment assigns", 'This edition assigns', 1),
               ('owner prose', 'editorial prose', 1)]}


def require_count(key: str, text: str, marker: str, expected: int = 1) -> None:
    actual = text.count(marker)
    if actual != expected:
        raise SystemExit(f"{key}: public transform expected {expected} occurrence(s) of {marker!r}, found {actual}")


def public_reading_copy(key: str, text: str) -> str:
    """Reduce private inputs only at reviewed, count-checked boundaries.

    This is intentionally not idempotent: cleaned snapshots are build inputs,
    never substitutes for private canonical inputs to --sync.
    """
    if key in {f"canon-{work_id.lower()}" for work_id, _, _ in DOSSIERS}:
        marker = "\n## Version receipt\n"
        require_count(key, text, marker)
        text = text.split(marker, 1)[0].rstrip() + "\n"
    for old, new, count in PUBLIC_REPLACEMENTS.get(key, []):
        require_count(key, text, old, count)
        text = text.replace(old, new)
    while "\n\n\n" in text:
        text = text.replace("\n\n\n", "\n\n")
    if key == "grounding":
        require_count(key, text, "---\n\n---")
        text = text.replace("---\n\n---", "---")
    return text


def public_payload(spec: dict, source_text: str) -> str:
    """Make only the documented public-reading reductions."""
    source_text = public_reading_copy(spec["key"], source_text)
    if spec.get("transform") == "protocol":
        marker = "\n## The fiber program, replayed\n"
        require_count(spec["key"], source_text, marker)
        return source_text.split(marker, 1)[0].rstrip() + "\n"
    if spec.get("transform") == "opening":
        for marker in ("<!-- GROUNDING:BEGIN ", "<!-- GROUNDING:END -->", "**Current work and its state:**"):
            require_count(spec["key"], source_text, marker)
        require_count(spec["key"], source_text, "<!-- GROUNDING:", 2)
        lines = [
            line
            for line in source_text.splitlines(keepends=True)
            if not line.startswith("<!-- GROUNDING:")
            and not line.startswith("**Current work and its state:**")
        ]
        return "".join(lines).strip() + "\n"
    return source_text


def sync_sources() -> None:
    SNAPSHOTS.mkdir(parents=True, exist_ok=True)
    rows = []
    recital = SNAPSHOTS / "00-recital.md"
    if not recital.is_file():
        raise SystemExit("canon/00-recital.md is missing")
    rows.append(
        {
            "key": "recital",
            "snapshot": recital.name,
            "canonical": "poams-exhibits-public/canon/00-recital.md",
            "sha256": sha(recital.read_bytes()),
            "bytes": recital.stat().st_size,
        }
    )
    pending = []
    for spec in SOURCES:
        target = SNAPSHOTS / spec["snapshot"]
        canonical = spec.get("canonical")
        canonical_sha = None
        if canonical:
            source = WORKSPACE / canonical
            if not source.is_file():
                raise SystemExit(f"missing canonical source: {source}")
            source_bytes = source.read_bytes()
            canonical_sha = sha(source_bytes)
            payload = public_payload(spec, source_bytes.decode()).encode()
            pending.append((target, payload))
        elif not target.is_file():
            raise SystemExit(f"missing maintained public source: {target}")
        if not canonical:
            payload = target.read_bytes()
        rows.append(
            {
                "key": spec["key"],
                "snapshot": spec["snapshot"],
                "canonical": canonical or f"poams-exhibits-public/canon/{spec['snapshot']}",
                "canonical_sha256": canonical_sha,
                "sha256": sha(payload),
                "bytes": len(payload),
            }
        )
    for target, payload in pending:
        target.write_bytes(payload)
    joined = "\n".join(f"{row['sha256']}  {row['snapshot']}" for row in rows) + "\n"
    manifest = {
        "schema": 2,
        "title": "The Canon source snapshot",
        "note": "canonical_sha256 identifies the workspace input; sha256 identifies the public payload",
        "sources": rows,
        "source_set_sha256": sha(joined.encode()),
    }
    (SNAPSHOTS / "SOURCE-MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    )
    (SNAPSHOTS / "SOURCE-MANIFEST.sha256").write_text(joined)


def load_manifest() -> dict:
    path = SNAPSHOTS / "SOURCE-MANIFEST.json"
    if not path.is_file():
        raise SystemExit("run scripts/build-canon-exhibit.py --sync first")
    manifest = json.loads(path.read_text())
    for row in manifest["sources"]:
        path = SNAPSHOTS / row["snapshot"]
        if not path.is_file() or sha(path.read_bytes()) != row["sha256"]:
            raise SystemExit(f"snapshot drift: {row['snapshot']}; rerun with --sync")
    return manifest


def render_markdown(path: Path, prefix: str) -> str:
    cmd = [
        "pandoc",
        str(path),
        "--from=gfm+tex_math_dollars",
        "--to=html5",
        "--mathjax",
        "--shift-heading-level-by=2",
        f"--id-prefix={prefix}-",
    ]
    return subprocess.run(cmd, check=True, text=True, capture_output=True).stdout.strip()


def source_card(key: str, status: str, payloads: dict[str, str]) -> str:
    source = next(row for row in load_manifest()["sources"] if row["key"] == key)
    path = SNAPSHOTS / source["snapshot"]
    payloads[key] = path.read_text()
    body = render_markdown(path, key)
    return f'''<article class="source-card" id="source-{html.escape(key)}">
  <div class="source-tools"><span>{html.escape(status)} · {source['bytes']:,} UTF-8 bytes</span><a href="canon/{source["snapshot"]}">Read Markdown</a><button type="button" data-copy-source="{html.escape(key)}" data-copy-label="{html.escape(key)}">Copy exact Markdown</button></div>
  <div class="source-document">{body}</div>
</article>'''


def group(number: str, anchor: str, title: str, lead: str, cards: list[str], opened: bool = False) -> str:
    statuses = {
        "core": "Compact Core · thirteen source-grounded dossiers · public edition",
        "grounding": "Grounding · adopted premises; Foundation Wall · amended October 18, 2026",
        "canon": "Source-grounded · thirteen-dossier public edition",
        "project-start-protocol": "Project-start gate · public excerpt",
        "opening-context-block": "Generated workspace context · unfilled editorial project fields",
        "canon-notes": "Editorial · public edition · scope and change control",
    }
    aliases = {"canon": "canon-pages", "project-start-protocol": "project-start", "opening-context-block": "opening-context"}
    alias = f'<span id="{aliases[anchor]}"></span>' if anchor in aliases else ""
    return f'''<details class="major" id="{anchor}"{' open' if opened else ''}>
  <summary><span class="major-number">{number}</span><div>{alias}<h2>{html.escape(title)}</h2><small>{html.escape(lead)}</small><small class="source-status">{html.escape(statuses[anchor])}</small></div></summary>
  <div class="major-body">{''.join(cards)}</div>
</details>'''


def listen_block() -> str:
    video = f'''<a class="pl-video" href="{html.escape(CANON_VIDEO_URL, quote=True)}" target="_blank" rel="noopener noreferrer" aria-label="Video in NotebookLM: The Canon (opens in a new tab)">VIDEO</a>'''
    binding_path = SNAPSHOTS / "audio-binding.json"
    if binding_path.is_file():
        binding = json.loads(binding_path.read_text())
        if binding.get("order") != 62 or not binding.get("audio_url"):
            raise SystemExit("canon/audio-binding.json has an invalid order or URL")
        body = f'''<div class="pl-pills"><button class="pl-pill" type="button" aria-label="Audio: The Canon" aria-expanded="false" aria-controls="pl-player-62">AUDIO</button>
{video}
</div>
<div id="pl-player-62" class="pl-player" hidden><button class="pl-play" type="button" aria-label="Play podcast" disabled>▶</button><input class="pl-seek" type="range" min="0" max="100" step="0.1" value="0" aria-label="Seek podcast" disabled><span class="pl-timer">0:00 / 0:00</span></div>
<p class="pl-error" role="status" hidden></p><audio preload="none" src="{html.escape(binding['audio_url'], quote=True)}"></audio>'''
    else:
        body = f'''<div class="pl-pills"><span class="pl-pill" aria-disabled="true" title="Audio under consideration">AUDIO</span>
{video}</div>'''
    wrapper = "poams-listen" if binding_path.is_file() else "canon-media-pending"
    return f'''<!-- POAMS LISTEN BEGIN -->
<div class="{wrapper}" data-podcast="62" data-podcast-title="The Canon">{body}</div>
<!-- POAMS LISTEN END -->'''


def build() -> None:
    manifest = load_manifest()
    payloads: dict[str, str] = {}
    recital_path = SNAPSHOTS / "00-recital.md"
    payloads["recital"] = recital_path.read_text()
    recital = render_markdown(recital_path, "recital")
    source_specs = {spec["key"]: spec for spec in SOURCES}
    canon_cards = [
        source_card(f"canon-{work_id.lower()}", source_specs[f"canon-{work_id.lower()}"]["status"], payloads)
        for work_id, _slug, _title in DOSSIERS
    ]
    sections = [
        f'''<section class="recital" id="recital"><div class="section-heading"><span>1</span><div><h2>Recital</h2><p>The reason this work comes first.</p></div></div><div class="source-tools"><span>Recital</span><a href="canon/00-recital.md">Read Markdown</a><button type="button" data-copy-source="recital" data-copy-label="recital">Copy exact Markdown</button></div><div class="source-document">{recital}</div></section>''',
        group("2", "core", "Core", "The compact thirteen-dossier orientation.", [source_card("canon-core", source_specs["canon-core"]["status"], payloads)], True),
        group("3", "grounding", "Grounding", "The program touchstone and its governing evidential wall.", [source_card("grounding", source_specs["grounding"]["status"], payloads), source_card("foundation-wall", source_specs["foundation-wall"]["status"], payloads)]),
        group("4", "canon", "Canon", "The authors’ words behind the thirteen source-grounded dossiers.", canon_cards),
        group("5", "project-start-protocol", "Project Start Protocol", "The gate before a construction begins.", [source_card("project-start-protocol", source_specs["project-start-protocol"]["status"], payloads)]),
        group("6", "opening-context-block", "Opening Context Block", "Legacy compatibility snapshot and project-specific fields; activation/migration against the thirteen-dossier edition is unreviewed.", [source_card("opening-context", source_specs["opening-context"]["status"], payloads), source_card("project-context", source_specs["project-context"]["status"], payloads)], True),
        group("7", "canon-notes", "Canon Notes", "Scope, status, open questions and change control.", [source_card("canon-notes", source_specs["canon-notes"]["status"], payloads)]),
    ]
    source_json = json.dumps(payloads, ensure_ascii=False).replace("<", "\\u003c")
    page = TEMPLATE
    replacements = {
        "__SOURCE_SET__": manifest["source_set_sha256"],
        "__LISTEN_BLOCK__": listen_block(),
        "__PLAYER_SCRIPT__": '<script defer src="./podcast-listen.js"></script>' if (SNAPSHOTS / "audio-binding.json").is_file() else "",
        "__SECTIONS__": "\n".join(sections),
        "__SOURCE_JSON__": source_json,
        "__COPY_ORDER__": json.dumps(["recital"] + [spec["key"] for spec in SOURCES]),
    }
    for old, new in replacements.items():
        page = page.replace(old, new)
    OUTPUT.write_text(page)


TEMPLATE = r'''<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<script id="canon-theme-init">
// Apply this page's preference before styles or body; never touch corpus theme keys.
(function(){let theme='dark';try{if(localStorage.getItem('canon-reading-theme')==='light')theme='light'}catch(e){}
document.documentElement.dataset.theme=theme;})();
</script>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta name="description" content="The source-grounded Normal Realism and POAMS canon, its core, grounding, project-start protocol, opening context and change notes.">
<title>The Canon | Normal Realism / POAMS</title>
<script>window.MathJax={tex:{inlineMath:[['$','$'],['\\(','\\)']],displayMath:[['$$','$$'],['\\[','\\]']]},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
<link rel="stylesheet" href="./reader-standard.css?v=20261011a">
<link rel="stylesheet" href="./podcast-listen.css">
__PLAYER_SCRIPT__
<style>
:root{--bg:#0d0d18;--card:#151523;--card2:#11111d;--edge:#2c2c42;--fg:#e8e4d8;--mut:#aaa28e;--gold:#d4af37;--blue:#5bc8f5;--green:#7bcf91;--shadow:0 18px 55px rgba(0,0,0,.28)}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;padding-top:72px!important;background:var(--bg);color:var(--fg);font-family:Georgia,serif;line-height:1.72;transition:background .25s,color .25s}main{max-width:1060px;margin:auto;padding:2rem 1.2rem 5rem;overflow-wrap:anywhere}a{color:var(--blue)}h1,h2,h3,h4{line-height:1.2}h1{color:var(--gold);font-size:clamp(2.25rem,7vw,4.6rem);font-weight:500;letter-spacing:-.035em;margin:.2rem 0 .45rem}h2{color:var(--gold);font-size:1.55rem}h3{color:var(--blue);font-size:1.22rem;margin-top:2rem}h4{color:var(--gold);font-size:1.05rem}strong,b{color:var(--blue)}em{color:var(--gold);font-style:normal;font-weight:600}code{font-size:.88em;background:rgba(127,127,127,.12);padding:.12em .32em;border-radius:4px}pre{overflow:auto;background:var(--card2);border:1px solid var(--edge);border-radius:8px;padding:1rem}blockquote{margin:1rem 0;padding:.15rem 0 .15rem 1rem;border-left:3px solid var(--gold);color:var(--mut)}
#poams-hub-nav{position:fixed;top:0;left:0;right:0;z-index:9999;height:44px;display:flex;align-items:center;justify-content:space-between;padding:0 1rem;background:rgba(13,13,24,.95);border-bottom:1px solid rgba(212,175,55,.25);font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;font-size:13px}#poams-hub-nav a{color:#d4af37;text-decoration:none;padding:5px 10px;white-space:nowrap}#poams-hub-nav .center{display:flex;flex-direction:column;align-items:center;line-height:1.2}#poams-hub-nav .title{font-weight:600;color:#e8e8f0;font-size:12px}#poams-hub-nav .pos{font-size:10px;color:rgba(212,175,55,.65);letter-spacing:.05em}
.hero{text-align:center;padding:2rem 0 1.2rem}.kicker{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;text-transform:uppercase;letter-spacing:.17em;color:var(--blue);font-size:.76rem}.subtitle{max-width:790px;margin:.7rem auto 1rem;color:var(--mut);font-size:1.05rem}.edition{font:600 .72rem/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:var(--mut);letter-spacing:.04em}.poams-listen,.canon-media-pending{margin:.9rem auto 0}.pl-video[aria-disabled="true"],.pl-pill[aria-disabled="true"]{opacity:.48;cursor:not-allowed}
.primer{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:.7rem;margin:1.5rem 0}.primer article{background:var(--card);border:1px solid var(--edge);border-radius:11px;padding:1rem}.primer h3{font:700 .82rem/1.3 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;text-transform:uppercase;letter-spacing:.08em;margin:0 0 .4rem}.primer p{font-size:.9rem;margin:0;color:var(--mut)}
.use-card{background:linear-gradient(135deg,rgba(212,175,55,.09),rgba(91,200,245,.05));border:1px solid rgba(212,175,55,.45);border-radius:14px;padding:1.2rem 1.3rem;margin:1.4rem 0}.use-card h3{margin:.05rem 0 .6rem}.use-card ol{margin:.4rem 0;padding-left:1.35rem}.copybar{display:flex;flex-wrap:wrap;gap:.55rem;margin:1rem 0 1.5rem}.copybar button,.source-tools button{border:1px solid var(--gold);background:transparent;color:var(--gold);border-radius:99px;padding:.48rem .8rem;font:650 .72rem/1.2 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;letter-spacing:.03em;cursor:pointer}.copybar button:hover,.source-tools button:hover{background:rgba(212,175,55,.1)}#copy-status{min-height:1.2rem;color:var(--green);font:600 .78rem/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.contents{display:flex;flex-wrap:wrap;gap:.45rem .6rem;margin:0 0 1.6rem}.contents a{border-bottom:1px solid rgba(91,200,245,.45);text-decoration:none;font:600 .76rem/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.recital,.major{background:var(--card);border:1px solid var(--edge);border-radius:14px;margin:1rem 0;box-shadow:var(--shadow)}.recital{padding:1.2rem 1.35rem}.section-heading{display:flex;align-items:center;gap:.8rem;border-bottom:1px solid var(--edge);padding-bottom:.8rem;margin-bottom:.7rem}.section-heading>span,.major-number{display:grid;place-items:center;min-width:2.25rem;height:2.25rem;border:1px solid var(--gold);border-radius:50%;color:var(--gold);font:700 .8rem/1 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.section-heading h2{margin:0}.section-heading p{margin:.15rem 0 0;color:var(--mut);font-size:.9rem}.major>summary{list-style:none;display:flex;align-items:center;gap:.8rem;padding:1rem 1.2rem;cursor:pointer}.major>summary::-webkit-details-marker{display:none}.major>summary::after{content:"+";margin-left:auto;color:var(--gold);font:400 1.5rem/1 sans-serif}.major[open]>summary::after{content:"−"}.major>summary h2{margin:0;display:block;color:var(--gold);font-size:1.2rem}.major>summary small{display:block;color:var(--mut);font-size:.84rem;margin-top:.15rem}.major-body{border-top:1px solid var(--edge);padding:0 1.25rem 1.25rem}.source-card{margin:1rem 0;padding:.2rem 0}.source-card+.source-card{border-top:1px dashed var(--edge);padding-top:1rem}.source-tools{display:flex;justify-content:space-between;align-items:center;gap:.7rem;color:var(--mut);font:600 .72rem/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.source-document{max-width:880px;margin:0 auto}.source-document>h3:first-child{color:var(--gold);font-size:1.4rem}.source-document li{margin:.28rem 0}.source-document table{display:block;overflow:auto;border-collapse:collapse}.source-document th,.source-document td{border:1px solid var(--edge);padding:.5rem;vertical-align:top}.source-document hr{border:0;border-top:1px solid var(--edge);margin:2rem 0}mjx-container[display="true"]{max-width:100%;overflow-x:auto;overflow-y:hidden}
.nav{display:flex;justify-content:space-between;gap:1rem;flex-wrap:wrap;font-size:.88rem;margin-top:2.5rem;border-top:1px solid var(--edge);padding-top:1rem}
@media(max-width:820px){.primer{grid-template-columns:repeat(2,minmax(0,1fr))}}@media(max-width:600px){#poams-hub-nav{padding-left:4px;padding-right:4px}#poams-hub-nav a{padding:5px 4px}.hide-sm{display:none}.primer{grid-template-columns:1fr}.source-tools{align-items:flex-start;flex-direction:column}.major>summary{padding:.9rem}.major-body,.recital{padding:.9rem}.source-document{font-size:.94rem}}
@media print{#poams-hub-nav,.poams-listen,.canon-media-pending,.copybar,#copy-status,.source-tools button{display:none!important}body{padding-top:0!important}details.major>summary::after{display:none}.major-body{border-top:1px solid #bbb}.recital,.major{box-shadow:none;break-inside:auto}}
[id]{scroll-margin-top:64px}.major>summary>div{min-width:0}.major-number{flex-shrink:0}.source-status{font-family:system-ui,sans-serif;font-size:.72rem!important}.canon-media-pending .pl-pills{display:flex;justify-content:center;gap:.6rem}.canon-media-pending .pl-pill,.canon-media-pending .pl-video{display:inline-block;border:1px solid var(--gold);border-radius:99px;padding:.3rem .8rem;color:var(--gold);font:600 .7rem/1.4 system-ui,sans-serif}.canon-media-pending a.pl-video{background:#6543a5;color:#fff;border-color:#6543a5;text-decoration:none}.canon-media-pending a.pl-video:hover{filter:brightness(1.15)}
/* Canon-only theme control. The standard dark palette and print rules stay intact. */
#canon-theme-toggle{min-height:44px;min-width:44px;border:1px solid var(--gold);border-radius:8px;background:var(--card);color:var(--gold);font:600 12px/1.2 system-ui,sans-serif;cursor:pointer;padding:.45rem .65rem;white-space:nowrap}
#canon-theme-toggle:focus-visible{outline:3px solid var(--blue);outline-offset:2px}
#canon-theme-toggle .theme-on{display:none}:root[data-theme="light"] #canon-theme-toggle .theme-off{display:none}:root[data-theme="light"] #canon-theme-toggle .theme-on{display:inline}
@media screen{
 html[data-theme="light"]:has(body.reader-standard),html[data-theme="light"] body.reader-standard{
  color-scheme:light!important;--bg:#f7f5ee!important;--card:#fffdf8!important;--card2:#eeeae0!important;--edge:#c9c3b4!important;
  --fg:#252735!important;--mut:#5d5f6b!important;--gold:#806115!important;--blue:#096184!important;--green:#27743d!important;--shadow:0 18px 55px rgba(45,38,20,.08);
 }
 html[data-theme="light"] body.reader-standard{background:var(--bg)!important;color:var(--fg)!important}
 html[data-theme="light"] body.reader-standard #poams-hub-nav{background:rgba(247,245,238,.96)!important;border-bottom-color:var(--edge)!important}
 html[data-theme="light"] body.reader-standard #poams-hub-nav a{color:var(--gold)!important}
 html[data-theme="light"] body.reader-standard #poams-hub-nav .title{color:var(--fg)!important}
 html[data-theme="light"] #poams-hub-nav .pos{color:var(--gold)}
}
@media print{#canon-theme-toggle{display:none!important}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
</style>
</head>
<body class="reader-standard">
<div id="poams-hub-nav"><a href="./index.html" aria-label="Previous: Hub">&larr; <span class="hide-sm">Hub</span></a><button type="button" id="canon-theme-toggle" aria-label="Light mode" aria-pressed="false">Light mode: <span class="theme-off">Off</span><span class="theme-on">On</span></button><div class="center"><a href="./index.html" style="display:flex;flex-direction:column;align-items:center;text-decoration:none"><span class="title">Return to Hub</span><span class="pos">NR · POAMS · CANON</span></a></div><a href="./the-rules-of-persistence.html" aria-label="Next: The Rules of Persistence"><span class="hide-sm">Rules of Persistence</span> &rarr;</a></div>
<main>
<header class="hero"><p class="kicker">Normal Realism first · POAMS follows</p><h1>The Canon</h1><p class="subtitle">Thirteen source-grounded dossiers, the program grounding built from them, and the compact tools that place both at the beginning of every new project.</p>__LISTEN_BLOCK__<p class="edition">Thirteen-dossier public edition &nbsp;·&nbsp; source set <code>__SOURCE_SET__</code></p></header>
<section class="primer" aria-label="Document roles"><article><h3>Canon</h3><p>The versioned, source-backed account of what Pope and Osborne hold.</p></article><article><h3>Grounding</h3><p>The canon applied with later adopted premises, provenance rules and measurement discipline.</p></article><article><h3>Opening context</h3><p>The compact generated package delivered into a working session. It does not create authority.</p></article><article><h3>Project start</h3><p>The gate that fixes the question, sources, premises, alternatives and failure test before construction.</p></article></section>
<section class="use-card"><h3>Workspace context</h3><p>This packet preserves workspace paths and <code>bin/corpus</code> commands. The bundled Markdown snapshots linked below are public reading copies of the frozen sources. Full source reading remains required; corpus locators are source locations, not public links. <a href="canon/SOURCE-MANIFEST.json">Source manifest</a>.</p><ol><li>The bundled <strong>Opening Context Block</strong> is a <strong>legacy compatibility snapshot</strong>. Its activation/migration against the thirteen-dossier edition is unreviewed; this public edition does not authorize activation.</li><li>For broad orientation, begin with the <strong>Canon Core</strong>.</li><li>For a project, use the order: <strong>Core → Grounding and Foundation Wall → relevant Canon dossier → Project-Start Protocol</strong>.</li><li>Use the full packet for reconciliation, amendment or review—not as the default context landslide.</li></ol><div class="copybar"><button type="button" id="copy-opening">Copy Workspace context</button><button type="button" data-copy-source="canon-core" data-copy-label="canon-core">Copy Canon Core</button><button type="button" id="copy-full">Copy Full Canon reference</button><button type="button" id="expand-all">Open all sections</button><button type="button" id="collapse-all">Collapse details</button></div><div id="copy-status" role="status" aria-live="polite"></div></section>
<nav class="contents" aria-label="Canon contents"><a href="#recital">1 · Recital</a><a href="#core">2 · Core</a><a href="#grounding">3 · Grounding</a><a href="#canon">4 · Canon</a><a href="#project-start-protocol">5 · Project Start Protocol</a><a href="#opening-context-block">6 · Opening Context Block</a><a href="#canon-notes">7 · Canon Notes</a></nav>
__SECTIONS__
<div class="nav"><a href="./index.html" aria-label="Previous: Hub">&larr; Hub</a><a href="./index.html">Return to the Hub</a><a href="./the-rules-of-persistence.html" aria-label="Next: The Rules of Persistence">The Rules of Persistence &rarr;</a></div>
</main>
<script id="canon-sources" type="application/json">__SOURCE_JSON__</script>
<script>
const themeButton=document.getElementById('canon-theme-toggle');
function updateThemeControl(){themeButton.setAttribute('aria-pressed',String(document.documentElement.dataset.theme==='light'));}
updateThemeControl();
themeButton.addEventListener('click',()=>{
  const theme=document.documentElement.dataset.theme==='dark'?'light':'dark';
  document.documentElement.dataset.theme=theme;
  try{localStorage.setItem('canon-reading-theme',theme)}catch(e){}
  updateThemeControl();
});
const canonSources=JSON.parse(document.getElementById('canon-sources').textContent);
const fullCanonOrder=__COPY_ORDER__;
async function copyText(text,label){
  let success=false;
  try{await navigator.clipboard.writeText(text);success=true}catch(e){
    const t=document.createElement('textarea');
    try{t.value=text;t.style.position='fixed';t.style.opacity='0';document.body.appendChild(t);t.select();success=document.execCommand('copy')===true}catch(e){success=false}finally{t.remove()}
  }
  const bytes=new TextEncoder().encode(text).length, characters=Array.from(text).length;
  document.getElementById('copy-status').textContent=success
    ? label+' copied: '+bytes+' UTF-8 bytes; '+characters+' Unicode characters.'
    : 'Unable to copy '+label+'. Use the bundled Markdown reading copies.';
  return success;
}
document.querySelectorAll('[data-copy-source]').forEach(b=>b.addEventListener('click',()=>copyText(canonSources[b.dataset.copySource],b.dataset.copyLabel)));
document.getElementById('copy-opening').addEventListener('click',()=>{const text=canonSources['opening-context'].trim()+'\n\n---\n\n'+canonSources['project-context'].trim()+'\n';copyText(text,'Workspace context')});
document.getElementById('copy-full').addEventListener('click',()=>{const text=fullCanonOrder.map(k=>'<!-- THE CANON · '+k.toUpperCase()+' -->\n\n'+canonSources[k].trim()).join('\n\n---\n\n')+'\n';copyText(text,'Full Canon reference')});
document.getElementById('expand-all').addEventListener('click',()=>document.querySelectorAll('details.major').forEach(d=>d.open=true));
document.getElementById('collapse-all').addEventListener('click',()=>document.querySelectorAll('details.major').forEach(d=>d.open=false));
let printStates=null;
addEventListener('beforeprint',()=>{
  if(printStates===null)printStates=Array.from(document.querySelectorAll('details.major'),d=>[d,d.open]);
  printStates.forEach(([d])=>d.open=true);
});
addEventListener('afterprint',()=>{
  if(printStates!==null){printStates.forEach(([d,open])=>d.open=open);printStates=null;}
});
let revealSequence=0;
async function revealTarget(){
  const sequence=++revealSequence;
  let id;try{id=decodeURIComponent(location.hash.slice(1))}catch(e){return}
  const target=id&&document.getElementById(id);if(!target)return;
  for(let node=target;node;node=node.parentElement)if(node.matches('details.major'))node.open=true;
  if(window.MathJax?.startup?.promise)await window.MathJax.startup.promise.catch(()=>{});
  await document.fonts.ready;
  await new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));
  if(sequence!==revealSequence)return;
  const header=document.getElementById('poams-hub-nav').getBoundingClientRect().bottom;
  window.scrollTo({top:Math.max(0,window.scrollY+target.getBoundingClientRect().top-header-16),behavior:matchMedia('(prefers-reduced-motion: reduce)').matches?'instant':'smooth'});
}
addEventListener('hashchange',revealTarget);addEventListener('load',revealTarget);revealTarget();
</script>
</body>
</html>
'''


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sync", action="store_true", help="refresh tracked snapshots from the workspace sources")
    args = parser.parse_args()
    if args.sync:
        sync_sources()
    build()
    print(f"built {OUTPUT} ({OUTPUT.stat().st_size} bytes)")
