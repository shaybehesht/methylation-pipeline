"""Shared MANGO branding: name, tagline, readability styling, and mango motifs."""
from __future__ import annotations

import streamlit as st

NAME = "MANGO"
TAGLINE = "Methylation Analysis for Novel Genomic Outcomes"
ICON = "🥭"

# Mango palette
_RIPE = "#F5A623"
_DEEP = "#E8590C"
_FLESH = "#FFD23F"
_LEAF = "#6DA34D"


def _slice_svg(filled: bool, size: int = 26) -> str:
    """A small mango-slice wedge SVG; ripe when filled, faint outline otherwise."""
    if filled:
        skin, flesh, stroke, opacity = _DEEP, _FLESH, _RIPE, "1"
    else:
        skin, flesh, stroke, opacity = "#EADFC7", "#F6EEDA", "#D9C79E", "0.7"
    return (
        f'<svg class="mango-slice" width="{size}" height="{size}" viewBox="0 0 40 40" '
        f'style="opacity:{opacity}" aria-hidden="true">'
        f'<path d="M6 30 C6 12 20 4 34 6 C30 20 20 32 6 30 Z" fill="{skin}"/>'
        f'<path d="M10 28 C11 16 20 9 30 10 C27 21 19 29 10 28 Z" fill="{flesh}"/>'
        f'<path d="M6 30 C6 12 20 4 34 6" fill="none" stroke="{stroke}" '
        f'stroke-width="1.6" stroke-linecap="round"/></svg>'
    )


_STYLE = f"""
<style>
/* Larger base size scales all rem-based text up for readability. */
html {{ font-size: 18px !important; }}
.stMarkdown p, .stMarkdown li {{ font-size: 1.07rem; line-height: 1.65; }}
h1 {{ font-size: 2.5rem !important; font-weight: 800 !important; letter-spacing: .2px; }}
h2 {{ font-size: 1.85rem !important; font-weight: 750 !important; }}
h3 {{ font-size: 1.45rem !important; font-weight: 700 !important; }}
/* Bold, larger widget labels (text inputs, selects, checkboxes, sliders, radios). */
[data-testid="stWidgetLabel"] p, label p {{ font-size: 1.08rem !important; font-weight: 650 !important; }}
/* Readable captions. */
[data-testid="stCaptionContainer"] p {{ font-size: 1.0rem !important; color: #6B4E2E !important; }}
/* Prominent buttons. */
.stButton button, .stDownloadButton button {{ font-size: 1.08rem !important; font-weight: 700 !important; }}
/* Bigger input/select/number text. */
.stTextInput input, .stNumberInput input, [data-baseweb="select"] {{ font-size: 1.05rem !important; }}
/* Metrics stand out. */
[data-testid="stMetricValue"] {{ font-size: 2rem !important; font-weight: 800 !important; }}
[data-testid="stMetricLabel"] p {{ font-size: 1.02rem !important; font-weight: 650 !important; }}
/* Sidebar navigation. */
[data-testid="stSidebarNav"] a span {{ font-size: 1.08rem !important; font-weight: 650 !important; }}
/* Alerts (info/warning/success) easier to read. */
[data-testid="stAlert"] p {{ font-size: 1.05rem !important; }}
/* Tables a touch larger. */
[data-testid="stDataFrame"], [data-testid="stTable"] {{ font-size: 1.02rem; }}

/* Mango accent bar + gentle per-tab entrance. */
.mango-bar {{ height: 8px; border-radius: 4px;
  background: linear-gradient(90deg, {_FLESH} 0%, {_RIPE} 45%, {_DEEP} 100%);
  margin-bottom: 0.85rem; animation: mango-rise .5s ease-out both; }}
h1 {{ animation: mango-rise .5s ease-out both; }}
@keyframes mango-rise {{ from {{ opacity: 0; transform: translateY(6px); }}
  to {{ opacity: 1; transform: translateY(0); }} }}

/* Slice meter. */
.mango-meter {{ display: flex; align-items: center; gap: 4px; margin: .2rem 0 .6rem; }}
.mango-meter .mango-slice {{ transition: opacity .3s ease; }}
.mango-meter .mango-slice-wrap.new .mango-slice {{ animation: mango-pop .45s ease-out both; }}
.mango-meter .cap {{ margin-left: .5rem; color: #6B4E2E; font-weight: 650; font-size: 1rem; }}
@keyframes mango-pop {{ 0% {{ transform: scale(.4); opacity: 0; }}
  60% {{ transform: scale(1.15); }} 100% {{ transform: scale(1); opacity: 1; }} }}

/* Running mango-slicing loader. */
.mango-loader {{ display: flex; align-items: center; gap: 6px; margin: .3rem 0 .2rem; }}
.mango-loader .mango-slice {{ animation: mango-slicing 1.6s ease-in-out infinite; }}
.mango-loader .mango-slice:nth-child(2) {{ animation-delay: .2s; }}
.mango-loader .mango-slice:nth-child(3) {{ animation-delay: .4s; }}
.mango-loader .mango-slice:nth-child(4) {{ animation-delay: .6s; }}
.mango-loader .lbl {{ margin-left: .5rem; color: {_DEEP}; font-weight: 700; }}
@keyframes mango-slicing {{ 0%, 100% {{ transform: translateY(0) rotate(0deg); opacity: .55; }}
  40% {{ transform: translateY(-6px) rotate(-8deg); opacity: 1; }} }}
</style>
<div class="mango-bar"></div>
"""


def style() -> None:
    """Inject MANGO readability styling, the mango accent bar, and animations."""
    st.markdown(_STYLE, unsafe_allow_html=True)


def slice_meter_html(completed: int, total: int, caption: str = "") -> str:
    """Return a row of mango slices: the first ``completed`` ripe, the rest faint."""
    total = max(0, int(total))
    completed = max(0, min(int(completed), total))
    slices = []
    for index in range(total):
        filled = index < completed
        # mark the most recently filled slice for the pop animation
        css_class = " new" if filled and index == completed - 1 else ""
        slices.append(f'<span class="mango-slice-wrap{css_class}">{_slice_svg(filled)}</span>')
    cap = f'<span class="cap">{caption}</span>' if caption else ""
    return f'<div class="mango-meter">{"".join(slices)}{cap}</div>'


def slice_meter(completed: int, total: int, caption: str = "") -> None:
    st.markdown(slice_meter_html(completed, total, caption), unsafe_allow_html=True)


def slicing_loader(label: str = "Slicing the mango…") -> str:
    """Return HTML for a looping mango-slicing loader (client-side CSS animation)."""
    wedges = "".join(_slice_svg(True, size=30) for _ in range(4))
    return f'<div class="mango-loader">{wedges}<span class="lbl">{label}</span></div>'


def accent_html(count: int = 3, size: int = 22) -> str:
    """Return a small decorative row of ripe mango slices."""
    slices = "".join(_slice_svg(True, size=size) for _ in range(max(0, int(count))))
    return f'<div class="mango-meter" style="margin-top:-0.5rem">{slices}</div>'


def accent(count: int = 3, size: int = 22) -> None:
    st.markdown(accent_html(count, size), unsafe_allow_html=True)
