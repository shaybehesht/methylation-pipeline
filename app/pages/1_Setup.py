from pathlib import Path

import streamlit as st

from app.file_picker import (
    BAM_EXTENSIONS,
    VCF_EXTENSIONS,
    data_roots,
    detect_bam_index,
    file_browser,
)
from app.state import initialize
from core.qc import gate, inspect_bam
from core.references import (
    available_assemblies,
    get_assembly,
    is_prepared,
    prepare_assembly,
    prepared_paths,
    reference_cache_root,
)

initialize()
st.title("🥭 1. Setup")
_roots = data_roots()
st.caption(
    "Select local files by browsing these locations: "
    + ", ".join(str(root) for root in _roots)
    + ". External drives are included automatically; nothing is uploaded through "
    "the browser and paths are never typed by hand."
)

samples = []
for index, existing in enumerate(st.session_state.samples):
    with st.expander(f"Sample {index + 1}", expanded=True):
        label = st.text_input(
            "Label", existing["label"], key=f"label_{index}",
            help="A unique display name used in comparisons and output filenames.",
        )
        bam = file_browser(
            "modBAM", key=f"bam_{index}", extensions=BAM_EXTENSIONS,
            help="Browse the mounted data root and select this sample's nanopore modBAM.",
        )
        if bam:
            index_path = detect_bam_index(bam)
            if index_path is not None:
                st.caption(f"Index detected: {Path(index_path).name}")
            else:
                st.warning(
                    "No BAM index (.bam.bai or .bai) was found next to this file. "
                    "Create one with `samtools index <file>.bam` so modkit can read "
                    "regions efficiently."
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
assembly_keys = available_assemblies()
current_assembly = st.session_state.get("assembly", assembly_keys[0])
if current_assembly not in assembly_keys:
    current_assembly = assembly_keys[0]
assembly = st.selectbox(
    "Genome build", assembly_keys,
    index=assembly_keys.index(current_assembly),
    format_func=lambda key: get_assembly(key).label,
    help="The managed FASTA, GENCODE annotation, and CpG islands for this build "
    "are downloaded once and cached for offline reuse.",
)
if assembly != st.session_state.get("assembly"):
    st.session_state.assembly = assembly
    st.session_state.reference_ready = False

st.caption(f"Reference cache: {reference_cache_root()}")

if is_prepared(assembly):
    paths = prepared_paths(assembly)
    st.session_state.reference_fasta = str(paths["fasta"])
    st.session_state.reference_gtf = str(paths["gtf"])
    st.session_state.reference_cpg_islands = str(paths["cpg_islands"])
    st.session_state.reference_ready = True
    st.success(f"{get_assembly(assembly).label} is cached and ready (offline).")
else:
    st.session_state.reference_ready = False
    st.info(
        f"{get_assembly(assembly).label} is not cached yet. Download and prepare it "
        "once; later runs are fully offline."
    )

if st.button(f"Download and prepare {get_assembly(assembly).label}"):
    progress_bar = st.progress(0.0)
    status = st.empty()

    def update(fraction: float, message: str) -> None:
        progress_bar.progress(min(max(fraction, 0.0), 1.0))
        status.write(message)

    try:
        paths = prepare_assembly(assembly, progress=update)
        st.session_state.reference_fasta = str(paths["fasta"])
        st.session_state.reference_gtf = str(paths["gtf"])
        st.session_state.reference_cpg_islands = str(paths["cpg_islands"])
        st.session_state.reference_ready = True
        st.success(f"{get_assembly(assembly).label} downloaded and prepared.")
    except Exception as exc:
        st.session_state.reference_ready = False
        st.error(f"Reference preparation failed: {exc}")

mod_options = ["5mC", "5hmC"]
current_mod = (st.session_state.get("modified_bases") or ["5mC"])[0]
if current_mod not in mod_options:
    current_mod = "5mC"
chosen_mod = st.selectbox(
    "Modified base", mod_options, index=mod_options.index(current_mod),
    help="Cytosine modification tabulated at CpG sites. modkit's --cpg pileup "
    "requires an explicit modified base; 5mC is standard for CpG methylation.",
)
st.session_state.modified_bases = [chosen_mod]
st.session_state.combine_strands = st.checkbox(
    "Combine CpG strands", value=bool(st.session_state.get("combine_strands", True)),
    help="Merge the two strands of each CpG (modkit --combine-strands). Matches "
    "the reference pipeline; disable only if the modBAM lacks MN tags.",
)

st.subheader("Optional phased VCF")
st.session_state.phased_vcf = file_browser(
    "Phased family VCF", key="phased_vcf", extensions=VCF_EXTENSIONS,
    help="Phased genotypes link alleles across variants and are needed to investigate "
    "methylation QTL and parent-of-origin effects. Leave unselected if unavailable.",
)
st.info("Optional fields improve interpretation. Leaving them blank does not prevent the basic trio DMR screen.")

if st.button("Validate BAMs and reference", type="primary"):
    try:
        if not st.session_state.get("reference_ready"):
            raise ValueError("Download and prepare the reference assembly first")
        if sum(item["role"] == "proband" for item in samples) != 1:
            raise ValueError("Select exactly one proband")
        if len({item["label"] for item in samples}) != 3:
            raise ValueError("Labels must be unique")
        for item in samples:
            if not item["bam_path"]:
                raise ValueError(f"Select a modBAM for {item['label']}")
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
