import json
from pathlib import Path

import pandas as pd
import streamlit as st

from app.state import initialize

initialize()
st.title("5. Results")
result = st.session_state.get("last_result")
if result is None:
    summary_path = Path(st.session_state.get("output_dir", "runs/latest")) / "summary.json"
    if summary_path.exists():
        result = json.loads(summary_path.read_text(encoding="utf-8"))

if result is None:
    st.info("Run an analysis to populate this page.")
    st.stop()

st.metric("Verdict", result["verdict"])
left, right = st.columns(2)
left.metric("Qualifying DMRs", result["candidate_count"])
right.metric("Candidate / null ratio", f"{result['ratio']:.2f}")
st.subheader("Interpretation")
st.write(result["reasoning"])

output = Path(result["output"])
figure = output / "dmr_effects.png"
if figure.exists():
    st.image(str(figure))
table_path = output / "proband_specific_DMRs.tsv"
if table_path.exists():
    candidates = pd.read_csv(table_path, sep="\t")
    st.subheader("Ranked candidates")
    st.dataframe(candidates, use_container_width=True)
    st.download_button(
        "Download TSV", candidates.to_csv(sep="\t", index=False),
        file_name="proband_specific_DMRs.tsv",
    )
report = Path(result["report"])
if report.exists():
    st.download_button(
        "Download self-contained HTML report", report.read_bytes(),
        file_name="methylation_trio_report.html", mime="text/html",
    )
