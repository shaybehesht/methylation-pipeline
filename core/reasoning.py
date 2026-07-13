"""Plain-language interpretation and self-contained report export."""
from __future__ import annotations

import base64
import html
from pathlib import Path

import pandas as pd


def explain(summary: dict, cutoff: float, caveats: list[str]) -> str:
    verdict = summary["verdict"]
    count, ratio = summary["candidate_count"], summary["ratio"]
    lead = {
        "PURSUE": "The screen supports follow-up.",
        "MARGINAL": "The screen is suggestive but not decisive.",
        "DO NOT PURSUE": "The screen does not currently support follow-up.",
    }[verdict]
    caveat_text = " ".join(caveats)
    return (
        f"{lead} {count} concordant proband-specific regions passed the empirical "
        f"null cutoff ({cutoff:.3g}); the candidate-to-null ratio was {ratio:.2f}. "
        f"This is an exploratory prioritization, not a diagnosis. {caveat_text}"
    )


def write_html_report(
    path: str | Path,
    title: str,
    summary: dict,
    reasoning: str,
    candidates: pd.DataFrame,
    figure: str | Path | None = None,
) -> Path:
    image = ""
    if figure and Path(figure).exists():
        encoded = base64.b64encode(Path(figure).read_bytes()).decode()
        image = f'<img alt="DMR effect plot" src="data:image/png;base64,{encoded}">'
    table = candidates.to_html(index=False, border=0, escape=True)
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem}}
table{{border-collapse:collapse;width:100%}}th,td{{padding:.4rem;border-bottom:1px solid #ddd}}
img{{max-width:100%}}.verdict{{font-size:1.5rem;font-weight:700}}</style></head>
<body><h1>{html.escape(title)}</h1><p class="verdict">{html.escape(str(summary["verdict"]))}</p>
<p>{html.escape(reasoning)}</p>{image}<h2>Ranked candidates</h2>{table}</body></html>"""
    output = Path(path)
    output.write_text(document, encoding="utf-8")
    return output
