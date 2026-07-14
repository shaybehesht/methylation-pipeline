import json
from pathlib import Path

import pandas as pd
import streamlit as st

from app.state import initialize
from core.reporting import build_run_archive

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
if result.get("evidence_status"):
    st.subheader("What can be assessed from the supplied metadata?")
    labels = {
        "phenotype": "Phenotype segregation",
        "parent_of_origin": "Parent of origin",
        "mqtl": "Methylation QTL",
        "tissue": "Tissue effects",
        "batch": "Batch effects",
    }
    evidence = pd.DataFrame([
        {"Question": labels.get(key, key), "Status": value.replace("_", " ")}
        for key, value in result["evidence_status"].items()
    ])
    st.dataframe(evidence, hide_index=True, use_container_width=True)
    st.caption(
        "“Inputs available” means the required metadata exists; it does not by "
        "itself establish causality or a diagnosis."
    )

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

st.subheader("Complete run archive")
st.caption(
    "The archive contains everything in the run directory: manifest, logs, "
    "bedMethyl pileups, pairwise DMR tables, ranked candidates, figures, and the "
    "HTML report. It can be large because pileups are included."
)
if st.button("Build complete ZIP"):
    with st.spinner("Building archive..."):
        archive = build_run_archive(output)
    size_mb = archive.stat().st_size / (1 << 20)
    st.session_state.archive_path = str(archive)
    st.success(f"Archive ready ({size_mb:.1f} MB): {archive}")
archive_path = st.session_state.get("archive_path")
if archive_path and Path(archive_path).exists() and Path(archive_path).parent == output:
    st.download_button(
        "Download complete ZIP", Path(archive_path).read_bytes(),
        file_name="methyl_trio_run.zip", mime="application/zip",
    )
