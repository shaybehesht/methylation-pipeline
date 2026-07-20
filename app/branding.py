"""Shared MANGO branding: name, tagline, and readability/theme styling."""
from __future__ import annotations

import streamlit as st

NAME = "MANGO"
TAGLINE = "Methylation Analysis for Novel Genomic Outcomes"
ICON = "🥭"

_STYLE = """
<style>
/* Larger base size scales all rem-based text up for readability. */
html { font-size: 18px !important; }
.stMarkdown p, .stMarkdown li { font-size: 1.07rem; line-height: 1.65; }
h1 { font-size: 2.5rem !important; font-weight: 800 !important; letter-spacing: .2px; }
h2 { font-size: 1.85rem !important; font-weight: 750 !important; }
h3 { font-size: 1.45rem !important; font-weight: 700 !important; }
/* Bold, larger widget labels (text inputs, selects, checkboxes, sliders, radios). */
[data-testid="stWidgetLabel"] p, label p { font-size: 1.08rem !important; font-weight: 650 !important; }
/* Readable captions. */
[data-testid="stCaptionContainer"] p { font-size: 1.0rem !important; color: #6B4E2E !important; }
/* Prominent buttons. */
.stButton button, .stDownloadButton button { font-size: 1.08rem !important; font-weight: 700 !important; }
/* Bigger input/select/number text. */
.stTextInput input, .stNumberInput input, [data-baseweb="select"] { font-size: 1.05rem !important; }
/* Metrics stand out. */
[data-testid="stMetricValue"] { font-size: 2rem !important; font-weight: 800 !important; }
[data-testid="stMetricLabel"] p { font-size: 1.02rem !important; font-weight: 650 !important; }
/* Sidebar navigation. */
[data-testid="stSidebarNav"] a span { font-size: 1.08rem !important; font-weight: 650 !important; }
/* Alerts (info/warning/success) easier to read. */
[data-testid="stAlert"] p { font-size: 1.05rem !important; }
/* Tables a touch larger. */
[data-testid="stDataFrame"], [data-testid="stTable"] { font-size: 1.02rem; }
/* Mango accent bar. */
.mango-bar { height: 8px; border-radius: 4px;
  background: linear-gradient(90deg, #FFD23F 0%, #F5A623 45%, #E8590C 100%);
  margin-bottom: 0.85rem; }
</style>
<div class="mango-bar"></div>
"""


def style() -> None:
    """Inject MANGO readability styling and the mango accent bar."""
    st.markdown(_STYLE, unsafe_allow_html=True)
