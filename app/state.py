"""Streamlit session-state helpers."""
from __future__ import annotations

import streamlit as st

from core.config import Affection, RegionConfig, Relationship, Role, Sample, Sex, TrioConfig
from core.thresholds import defaults


def initialize() -> None:
    st.session_state.setdefault("samples", [
        {"label": "Proband", "bam_path": "", "sex": "F", "role": "proband", "affection": "affected"},
        {"label": "Relative 1", "bam_path": "", "sex": "F", "role": "relative"},
        {"label": "Relative 2", "bam_path": "", "sex": "M", "role": "relative"},
    ])
    st.session_state.setdefault("assembly", "hg38")
    st.session_state.setdefault("modified_bases", ["5mC"])
    st.session_state.setdefault("combine_strands", True)
    st.session_state.setdefault("reference_fasta", "")
    st.session_state.setdefault("reference_gtf", "")
    st.session_state.setdefault("reference_cpg_islands", "")
    st.session_state.setdefault("reference_ready", False)
    st.session_state.setdefault("region_mode", "chromosomes")
    st.session_state.setdefault("chromosomes", ["chr1", "chr2", "chr11", "chr15"])
    st.session_state.setdefault("genes", ["MECP2", "UBE3A", "SNRPN"])
    st.session_state.setdefault("thresholds", defaults())
    st.session_state.setdefault("qc_passed", False)
    st.session_state.setdefault("phased_vcf", "")


def config() -> TrioConfig:
    initialize()
    samples = [
        Sample(
            item["label"], item["bam_path"], Sex(item["sex"]), Role(item["role"]),
            Relationship(item["relationship"]) if item.get("relationship") else None,
            Affection(item["affection"]) if item.get("affection") else None,
        )
        for item in st.session_state.samples
    ]
    return TrioConfig(
        samples=samples,
        reference_fasta=st.session_state.reference_fasta,
        output_dir=st.session_state.get("output_dir", "runs/latest"),
        regions=RegionConfig(
            mode=st.session_state.region_mode,
            chromosomes=st.session_state.chromosomes,
            genes=st.session_state.genes,
        ),
        thresholds=st.session_state.thresholds,
        phased_vcf=st.session_state.phased_vcf,
        assembly=st.session_state.assembly,
        modified_bases=list(st.session_state.modified_bases),
        combine_strands=bool(st.session_state.get("combine_strands", True)),
    )
