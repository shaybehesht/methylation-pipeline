"""Methylation Trio Platform entry point."""
import streamlit as st

from app.state import initialize

st.set_page_config(page_title="Methylation Trio", page_icon="🧬", layout="wide")
initialize()

st.title("Methylation Trio Platform")
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
    "1. **Setup** — browse local BAMs, roles, sexes, and a managed reference assembly; run QC.\n"
    "2. **Regions** — choose whole genome, chromosomes, or target genes.\n"
    "3. **Thresholds** — review and edit every analysis cutoff.\n"
    "4. **Run** — execute modkit locally with progress and logs.\n"
    "5. **Results** — review candidates, reasoning, plots, and the HTML report."
)
