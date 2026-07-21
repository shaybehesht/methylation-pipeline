"""MANGO — Methylation Analysis for Novel Genomic Outcomes: entry point."""
import streamlit as st

from app import branding
from app.state import initialize

st.set_page_config(page_title=branding.NAME, page_icon=branding.ICON, layout="wide")
initialize()
branding.style()

st.title(f"{branding.ICON} {branding.NAME}")
st.caption(branding.TAGLINE)
branding.accent(3)
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
    "0. **AnVIL data** *(optional)* — pull modBAMs from your AnVIL / Terra "
    "workspace bucket on Google Cloud, or **Remote data** to mount BAMs from an "
    "SSH server.\n"
    "1. **Setup** — browse local BAMs, roles, sexes, and a managed reference assembly; run QC.\n"
    "2. **Regions** — choose whole genome, chromosomes, or target genes.\n"
    "3. **Thresholds** — review and edit every analysis cutoff.\n"
    "4. **Run** — execute modkit locally with progress and logs.\n"
    "5. **Results** — review candidates, reasoning, plots, and the HTML report."
)
