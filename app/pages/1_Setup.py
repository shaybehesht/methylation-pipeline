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
        label = st.text_input(
            "Label", existing["label"], key=f"label_{index}",
            help="A unique display name used in comparisons and output filenames.",
        )
        bam = st.text_input(
            "BAM path", existing["bam_path"], key=f"bam_{index}", placeholder="/data/sample.bam",
            help="Path to this sample's nanopore modBAM. Required to run analysis.",
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
        tissue = st.text_input(
            "Tissue type (optional)", existing.get("tissue", ""), key=f"tissue_{index}",
            placeholder="e.g. blood",
            help="DNA methylation is tissue-specific. Different or missing tissues limit interpretation.",
        )
        batch = st.text_input(
            "Sequencing/library batch (optional)", existing.get("batch", ""), key=f"batch_{index}",
            placeholder="e.g. flowcell or library ID",
            help="Differences in library preparation, flowcell, run, or basecalling can mimic biological methylation differences.",
        )
        samples.append({
            "label": label, "bam_path": bam, "sex": sex, "role": role,
            "relationship": None if relationship == "Not provided" else relationship,
            "affection": None if affection == "Not provided" else affection,
            "tissue": tissue, "batch": batch,
        })
st.session_state.samples = samples
st.session_state.reference_fasta = st.text_input(
    "Reference FASTA", st.session_state.reference_fasta, placeholder="/data/reference.fa",
    help="The exact genome assembly used to align all BAMs. Contig names and lengths are checked.",
)
st.session_state.phased_vcf = st.text_input(
    "Phased family VCF (optional)", st.session_state.phased_vcf,
    placeholder="/data/family.phased.vcf.gz",
    help="Phased genotypes link alleles across variants and are needed to investigate methylation QTL and parent-of-origin effects. Leave blank if unavailable.",
)
st.info("Optional fields improve interpretation. Leaving them blank does not prevent the basic trio DMR screen.")

if st.button("Validate BAMs and reference", type="primary"):
    try:
        if sum(item["role"] == "proband" for item in samples) != 1:
            raise ValueError("Select exactly one proband")
        if len({item["label"] for item in samples}) != 3:
            raise ValueError("Labels must be unique")
        if not Path(st.session_state.reference_fasta).exists():
            raise FileNotFoundError(st.session_state.reference_fasta)
        if st.session_state.phased_vcf and not Path(st.session_state.phased_vcf).exists():
            raise FileNotFoundError(st.session_state.phased_vcf)
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
