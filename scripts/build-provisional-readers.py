#!/usr/bin/env python3
"""Build media-header readers around the unchanged public provisional PDFs."""
import json
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ROWS = json.loads((ROOT / "provisionals/media.json").read_text())["provisionals"]


def audio_control(row):
    title = escape(f'{row["id"]}: {row["title"]}', quote=True)
    if not row["audio"]:
        return ('<button class="pl-pill pl-pending" type="button" disabled '
                'aria-label="Audio not yet available">AUDIO <small>Coming soon</small></button>'), ""
    audio = row["audio"]
    if not audio.startswith("https://audio.resonantmachines.us/"):
        raise ValueError(f'{row["id"]}: unapproved audio host')
    ident = row["id"].lower()
    control = (f'<button class="pl-pill" type="button" aria-label="Audio: {title}" '
               f'aria-expanded="false" aria-controls="pl-player-{ident}">AUDIO</button>')
    panel = f'''<div id="pl-player-{ident}" class="pl-player" hidden><button class="pl-play" type="button" aria-label="Play podcast" disabled>▶</button><input class="pl-seek" type="range" min="0" max="100" step="0.1" value="0" aria-label="Seek podcast" disabled><span class="pl-timer">0:00 / 0:00</span></div>
<p class="pl-error" role="status" hidden></p><audio preload="none" src="{escape(audio, quote=True)}"></audio>'''
    return control, panel


for index, row in enumerate(ROWS):
    ident = row["id"]
    if not (ROOT / "provisionals" / f"{ident}.pdf").is_file():
        raise FileNotFoundError(f"Missing original {ident}.pdf")
    title = escape(row["title"])
    label = escape(row["document_label"])
    video = escape(row["video"], quote=True)
    audio_pill, audio_panel = audio_control(row)
    playback = '<script defer src="../podcast-listen.js"></script>' if row["audio"] else ""
    previous = (f'<a rel="prev" href="{ROWS[index-1]["id"]}.html">← {ROWS[index-1]["id"]}</a>'
                if index else '<a href="../index.html#patents-pending">← All Provisionals</a>')
    following = (f'<a rel="next" href="{ROWS[index+1]["id"]}.html">{ROWS[index+1]["id"]} →</a>'
                 if index + 1 < len(ROWS) else '<a href="../index.html#patents-pending">All Provisionals →</a>')
    html = f'''<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{ident} · {title} — Normal Realism &amp; POAMS</title>
<meta name="description" content="{ident}: {title}. Read the provisional application and explore its companion media.">
<link rel="stylesheet" href="../podcast-listen.css">
<link rel="stylesheet" href="reader.css">
{playback}
</head>
<body>
<a class="skip-link" href="#document">Skip to document</a>
<nav class="reader-nav" aria-label="Resource Hub"><a href="../index.html#patents-pending">← Return to Hub</a><span class="family">Normal Realism · POAMS</span></nav>
<header class="reader-header">
<p class="eyebrow">Provisional Application · {ident}</p>
<h1>{title}</h1>
<p class="document-meta">{label} · {row["pages"]} pages</p>
<div class="poams-listen" data-provisional="{ident}">
<div class="pl-pills">{audio_pill}<a class="pl-video" href="{video}" target="_blank" rel="noopener noreferrer" aria-label="Video in NotebookLM: {ident}, {title} (opens in a new tab)">VIDEO</a></div>
{audio_panel}
</div>
<div class="document-actions"><a href="{ident}.pdf" target="_blank" rel="noopener" type="application/pdf">Open PDF ↗</a><a href="{ident}.pdf" download>Download PDF ↓</a></div>
</header>
<main id="document" class="document-shell" tabindex="-1">
<iframe class="pdf-reader" src="{ident}.pdf#view=FitH" title="{ident}: {title} — {row["pages"]}-page PDF" referrerpolicy="no-referrer"></iframe>
<p class="reader-help">Prefer your PDF reader? <a href="{ident}.pdf" target="_blank" rel="noopener">Open the document in a new tab.</a></p>
</main>
<nav class="reader-pagination" aria-label="Provisional applications">{previous}<span class="position">{index+1} / {len(ROWS)}</span>{following}</nav>
</body>
</html>
'''
    (ROOT / "provisionals" / f"{ident}.html").write_text(html)

print(f"Built {len(ROWS)} provisional readers; original PDFs untouched.")
