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
    if summary.get("design") == "phenotype_segregation":
        finding = (
            f"{count} regions were similar between affected samples and differed "
            f"concordantly from the unaffected relative at the {cutoff:.3g} effect cutoff"
        )
    else:
        finding = (
            f"{count} concordant proband-specific regions passed the empirical "
            f"null cutoff ({cutoff:.3g})"
        )
    denominator = summary.get("denominator_label", "null regions")
    return (
        f"{lead} {finding}; the candidate-to-{denominator} ratio was {ratio:.2f}. "
        f"This is an exploratory prioritization, not a diagnosis. {caveat_text}"
    )


def write_html_report(
    path: str | Path,
    title: str,
    summary: dict,
    reasoning: str,
    candidates: pd.DataFrame,
    figure: str | Path | None = None,
    figures: list[str | Path] | None = None,
) -> Path:
    sources = figures if figures is not None else ([figure] if figure else [])
    images = []
    for source in sources:
        if source and Path(source).exists():
            encoded = base64.b64encode(Path(source).read_bytes()).decode()
            caption = html.escape(Path(source).stem.replace("_", " "))
            images.append(
                f'<figure><img alt="{caption}" src="data:image/png;base64,{encoded}">'
                f'<figcaption>{caption}</figcaption></figure>'
            )
    image = "\n".join(images)
    table = candidates.to_html(index=False, border=0, escape=True)
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#3D2B1F}}
.mango-bar{{height:8px;border-radius:4px;background:linear-gradient(90deg,#FFD23F 0%,#F5A623 45%,#E8590C 100%);margin-bottom:1rem}}
h1,h2{{color:#B45309}}
table{{border-collapse:collapse;width:100%}}th{{background:#FFEFC9;text-align:left}}
th,td{{padding:.4rem;border-bottom:1px solid #F3E0B5}}
img{{max-width:100%}}.verdict{{font-size:1.5rem;font-weight:700;color:#E8590C}}
figure{{margin:1.2rem 0}}figcaption{{color:#8A6D3B;font-size:.9rem}}</style></head>
<body><div class="mango-bar"></div>
<h1>🥭 {html.escape(title)}</h1><p class="verdict">{html.escape(str(summary["verdict"]))}</p>
<p>{html.escape(reasoning)}</p>{image}<h2>Ranked candidates</h2>{table}</body></html>"""
    output = Path(path)
    output.write_text(document, encoding="utf-8")
    return output
