from pathlib import Path

import streamlit as st

from app.state import initialize
from core.qc import gate, inspect_bam

initialize()
st.title("1. Setup")
st.caption("Paths are resolved inside the container; mount host data at /data.")

samples = []
for index, existing in enumerate(st.session_state.samples):
    with st.expander(f"Sample {index + 1}", expanded=True):
        label = st.text_input("Label", existing["label"], key=f"label_{index}")
        bam = st.text_input("BAM path", existing["bam_path"], key=f"bam_{index}", placeholder="/data/sample.bam")
        sex = st.radio("Sex", ["F", "M"], index=["F", "M"].index(existing["sex"]), horizontal=True, key=f"sex_{index}")
        role = st.selectbox("Role", ["proband", "relative"], index=["proband", "relative"].index(existing["role"]), key=f"role_{index}")
        samples.append({"label": label, "bam_path": bam, "sex": sex, "role": role})
st.session_state.samples = samples
st.session_state.reference_fasta = st.text_input(
    "Reference FASTA", st.session_state.reference_fasta, placeholder="/data/reference.fa"
)

if st.button("Validate BAMs and reference", type="primary"):
    try:
        if sum(item["role"] == "proband" for item in samples) != 1:
            raise ValueError("Select exactly one proband")
        if len({item["label"] for item in samples}) != 3:
            raise ValueError("Labels must be unique")
        if not Path(st.session_state.reference_fasta).exists():
            raise FileNotFoundError(st.session_state.reference_fasta)
        results = [
            inspect_bam(item["bam_path"], st.session_state.reference_fasta)
            for item in samples
        ]
        passed, errors = gate(results)
        st.session_state.qc_results = results
        st.session_state.qc_passed = passed
        for result in results:
            st.write(
                f"**{Path(result['bam']).name}** — model: {result['basecaller_model']}; "
                f"HP-tagged reads: {result['hp_fraction']:.1%}; reference match: "
                f"{'yes' if result['reference_matches'] else 'no'}"
            )
            if not result["has_hp_tags"]:
                st.warning("No HP tags detected. Haplotype read plots will be unavailable.")
        if errors:
            for error in errors:
                st.error(error)
        else:
            st.success("Required QC gates passed.")
    except Exception as exc:
        st.session_state.qc_passed = False
        st.error(str(exc))
