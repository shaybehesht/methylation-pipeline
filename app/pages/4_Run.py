from datetime import datetime, timezone
from pathlib import Path

import streamlit as st

from app.state import config, initialize
from core.pipeline import run

initialize()
st.title("4. Run")
default_output = f"runs/run-{datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')}"
st.session_state.output_dir = st.text_input(
    "Output directory", st.session_state.get("output_dir", default_output)
)

gtf = st.session_state.get("reference_gtf") or None
cpg_islands = st.session_state.get("reference_cpg_islands") or None
if gtf:
    st.caption(f"GENCODE annotation (from selected assembly): {gtf}")

if not st.session_state.qc_passed:
    st.warning("Run Setup QC successfully before analysis.")

if st.button("Run methylation analysis", type="primary", disabled=not st.session_state.qc_passed):
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
        st.session_state.last_result = result
        st.success(f"Complete: {result['verdict']}")
        log = Path(result["output"]) / "pipeline.log"
        if log.exists():
            with st.expander("Pipeline log"):
                st.code(log.read_text(encoding="utf-8", errors="replace"))
    except Exception as exc:
        st.error(f"Analysis failed: {exc}")
        output = Path(st.session_state.output_dir)
        log = output / "pipeline.log"
        if log.exists():
            with st.expander("Pipeline log", expanded=True):
                st.code(log.read_text(encoding="utf-8", errors="replace"))
