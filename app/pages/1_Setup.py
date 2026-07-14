from pathlib import Path

import streamlit as st

from app.file_picker import bam_index, pick_file
from app.state import initialize
from core.qc import gate, inspect_bam
from core.references import ASSEMBLIES, ensure_assembly, is_ready

initialize()
st.title("1. Setup")
st.caption(
    "Select files from the mounted data directory. Nothing is uploaded through "
    "the browser; the app reads files in place."
)

samples = []
for index, existing in enumerate(st.session_state.samples):
    with st.expander(f"Sample {index + 1}", expanded=True):
        label = st.text_input(
            "Label", existing["label"], key=f"label_{index}",
            help="A unique display name used in comparisons and output filenames.",
        )
        bam = pick_file(
            "BAM file", key=f"bam_{index}", extensions=(".bam",),
            help="This sample's nanopore modBAM, chosen from the data directory.",
        )
        bam = bam or existing.get("bam_path", "")
        if bam:
            if bam_index(bam):
                st.caption(f"Index found: {Path(bam_index(bam)).name}")
            else:
                st.warning(
                    "No .bai index found next to this BAM. Create one with "
                    "`samtools index` so region fetches work."
                )
        sex = st.radio(
            "Sex", ["F", "M"], index=["F", "M"].index(existing["sex"]),
            horizontal=True, key=f"sex_{index}",
            help="Used only to decide whether chrX or chrY comparisons are technically valid.",
        )
        role = st.selectbox(
            "Analysis role", ["proband", "relative"],
            index=["proband", "relative"].index(existing["role"]), key=f"role_{index}",
            help="Exactly one sample must be the proband; the other two are family comparators.",
        )
        relationship_options = ["Not provided", "mother", "father", "sibling", "other"]
        relationship_value = existing.get("relationship") or "Not provided"
        relationship = st.selectbox(
            "Relationship (optional)", relationship_options,
            index=relationship_options.index(relationship_value), key=f"relationship_{index}",
            disabled=role == "proband",
            help="Identifies parental lineage. Both mother and father plus a phased VCF are needed to assign parent of origin.",
        )
        affection_options = ["Not provided", "affected", "unaffected", "unknown"]
        affection_value = existing.get("affection") or "Not provided"
        affection = st.selectbox(
            "Clinical status (optional)", affection_options,
            index=affection_options.index(affection_value), key=f"affection_{index}",
            help="Affected samples are expected to resemble one another at disease-associated DMRs; unaffected samples provide contrast. Unknown is recorded uncertainty.",
        )
        samples.append({
            "label": label, "bam_path": bam, "sex": sex, "role": role,
            "relationship": None if relationship == "Not provided" else relationship,
            "affection": None if affection == "Not provided" else affection,
        })
st.session_state.samples = samples

st.subheader("Reference assembly")
assembly_keys = list(ASSEMBLIES)
st.session_state.assembly = st.selectbox(
    "Genome assembly", assembly_keys,
    index=assembly_keys.index(st.session_state.assembly),
    format_func=lambda key: ASSEMBLIES[key].label,
    help="The assembly your BAMs were aligned to. It is downloaded and cached once, then reused offline.",
)
if is_ready(st.session_state.assembly):
    st.success(f"{ASSEMBLIES[st.session_state.assembly].label} is prepared and cached.")
else:
    st.info("This assembly is not cached yet. Download and prepare it once before running.")
if st.button("Download and prepare reference"):
    progress_bar = st.progress(0.0)
    status = st.empty()

    def update(fraction: float, message: str) -> None:
        progress_bar.progress(min(max(fraction, 0.0), 1.0))
        status.write(message)

    try:
        bundle = ensure_assembly(st.session_state.assembly, progress=update)
        st.session_state.reference_fasta = bundle.fasta
        st.session_state.reference_gtf = bundle.gtf
        st.session_state.reference_cpg = bundle.cpg_islands
        st.session_state.reference_ready = True
        st.success("Reference prepared and cached.")
    except Exception as exc:
        st.session_state.reference_ready = False
        st.error(f"Reference preparation failed: {exc}")

st.subheader("Optional phased VCF")
vcf = pick_file(
    "Phased family VCF", key="phased_vcf", extensions=(".vcf.gz", ".vcf", ".bcf"),
    help="Phased genotypes are needed to investigate methylation QTL and parent-of-origin effects. Optional.",
)
st.session_state.phased_vcf = vcf or ""

st.info("Optional fields improve interpretation. Leaving them blank does not prevent the basic trio DMR screen.")

if st.button("Validate BAMs and reference", type="primary"):
    try:
        if sum(item["role"] == "proband" for item in samples) != 1:
            raise ValueError("Select exactly one proband")
        if len({item["label"] for item in samples}) != 3:
            raise ValueError("Labels must be unique")
        if any(not item["bam_path"] for item in samples):
            raise ValueError("Select a BAM file for every sample")
        if not is_ready(st.session_state.assembly):
            raise ValueError("Download and prepare the reference assembly first")
        reference = st.session_state.reference_fasta
        results = [inspect_bam(item["bam_path"], reference) for item in samples]
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
