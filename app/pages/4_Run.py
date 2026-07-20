from datetime import datetime, timezone
from pathlib import Path

import streamlit as st

from app import branding
from app.file_picker import data_roots
from app.state import config, initialize

from core.pipeline import run

initialize()
branding.style()
st.title("🥭 4. Run")

st.subheader("Where to save this run")
run_name = st.text_input(
    "Run name", st.session_state.get("run_name", f"run-{datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')}"),
    help="A subfolder created inside the chosen location.",
)
st.session_state.run_name = run_name
location_options = ["Project ./runs"] + [str(root) for root in data_roots()]
chosen_location = st.selectbox(
    "Save location", location_options,
    help="Write results to the project's ./runs folder or to any browsable data "
    "root, including an external drive.",
)
if chosen_location == "Project ./runs":
    base = Path("runs")
else:
    base = Path(chosen_location) / "methyl_trio_runs"
st.session_state.output_dir = str(base / run_name)
st.caption(f"Output directory: {st.session_state.output_dir}")

mode_labels = {
    "whole_genome": "Genome-wide HMM segmentation (all autosomes + chrX)",
    "chromosomes": "Per-chromosome HMM segmentation (selected chromosomes)",
    "targeted": "Targeted per-region paired Wilcoxon (gene panel)",
}
st.caption(f"Analysis: {mode_labels.get(st.session_state.get('region_mode'), st.session_state.get('region_mode'))}")

gtf = st.session_state.get("reference_gtf") or None
cpg_islands = st.session_state.get("reference_cpg_islands") or None
if gtf:
    st.caption(f"GENCODE annotation (from selected assembly): {gtf}")

if not st.session_state.qc_passed:
    st.warning("Run Setup QC successfully before analysis.")

if st.button("Run methylation analysis", type="primary", disabled=not st.session_state.qc_passed):
    loader = st.empty()
    loader.markdown(branding.slicing_loader(), unsafe_allow_html=True)
    progress_bar = st.progress(0.0)
    status = st.empty()

    def update(fraction: float, message: str) -> None:
        progress_bar.progress(fraction)
        status.write(message)

    try:
        run_kwargs = {"gtf": gtf, "progress": update}
        if cpg_islands:
            run_kwargs["cpg_islands"] = cpg_islands
        result = run(config(), **run_kwargs)
        loader.empty()
        st.session_state.last_result = result
        st.success(f"Complete: {result['verdict']}")
        log = Path(result["output"]) / "pipeline.log"
        if log.exists():
            with st.expander("Pipeline log"):
                st.code(log.read_text(encoding="utf-8", errors="replace"))
    except Exception as exc:
        loader.empty()
        st.error(f"Analysis failed: {exc}")
        output = Path(st.session_state.output_dir)
        log = output / "pipeline.log"
        if log.exists():
            with st.expander("Pipeline log", expanded=True):
                st.code(log.read_text(encoding="utf-8", errors="replace"))
