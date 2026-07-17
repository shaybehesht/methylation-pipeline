"""MANGO — Methylation Analysis for Novel Genomic Outcomes: entry point."""
import streamlit as st

from app.state import initialize

st.set_page_config(page_title="MANGO", page_icon="🥭", layout="wide")
initialize()

st.markdown(
    """
    <style>
    .mango-bar {
        height: 8px;
        background: linear-gradient(90deg, #FFD23F 0%, #F5A623 45%, #E8590C 100%);
        border-radius: 4px;
        margin-bottom: 0.75rem;
    }
    </style>
    <div class="mango-bar"></div>
    """,
    unsafe_allow_html=True,
)

st.title("🥭 MANGO")
st.caption("Methylation Analysis for Novel Genomic Outcomes")
st.write(
    "Offline, configurable three-sample nanopore methylation analysis. "
    "Use the pages in the sidebar to validate inputs, define the region scope, "
    "adjust documented thresholds, run the analysis, and inspect results."
)
st.info(
    "This is an exploratory family screen. It cannot establish diagnosis or "
    "population-level significance from three samples."
)
st.subheader("Workflow")
st.markdown(
    "0. **Remote data** *(optional)* — mount BAMs from an SSH server to browse them by path.\n"
    "1. **Setup** — browse local BAMs, roles, sexes, and a managed reference assembly; run QC.\n"
    "2. **Regions** — choose whole genome, chromosomes, or target genes.\n"
    "3. **Thresholds** — review and edit every analysis cutoff.\n"
    "4. **Run** — execute modkit locally with progress and logs.\n"
    "5. **Results** — review candidates, reasoning, plots, and the HTML report."
)
