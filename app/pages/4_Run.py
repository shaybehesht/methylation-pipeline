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
gtf = st.text_input("GENCODE GTF", "/app/annotations/gencode.annotation.gtf.gz")

if not st.session_state.qc_passed:
    st.warning("Run Setup QC successfully before analysis.")

if st.button("Run methylation analysis", type="primary", disabled=not st.session_state.qc_passed):
    progress_bar = st.progress(0.0)
    status = st.empty()

    def update(fraction: float, message: str) -> None:
        progress_bar.progress(fraction)
        status.write(message)

    try:
        result = run(config(), gtf=gtf, progress=update)
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
