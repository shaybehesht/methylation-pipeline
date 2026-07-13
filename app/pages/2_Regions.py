import streamlit as st

from app.state import initialize

initialize()
st.title("2. Regions")
labels = {
    "whole_genome": "Whole genome",
    "chromosomes": "Selected chromosomes",
    "targeted": "Targeted gene panel",
}
keys = list(labels)
mode = st.radio(
    "Analysis scope", keys, format_func=labels.get,
    index=keys.index(st.session_state.region_mode),
)
st.session_state.region_mode = mode

if mode == "chromosomes":
    autosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    st.session_state.chromosomes = st.multiselect(
        "Chromosomes", autosomes, default=st.session_state.chromosomes
    )
    st.caption("chrX is retained only for female–female comparisons; chrY only for male–male. chrM is always excluded.")
elif mode == "targeted":
    panel_a = "MECP2, UBE3A, SNRPN, MAGEL2, SNORD116"
    panel_b = "FMR1, SHANK3, CHD8, EHMT1, KMT2D"
    preset = st.selectbox("Starting panel", ["Current", "Panel A", "Panel B"])
    initial = ", ".join(st.session_state.genes)
    if preset == "Panel A":
        initial = panel_a
    elif preset == "Panel B":
        initial = panel_b
    pasted = st.text_area("Gene symbols (comma, space, or newline separated)", initial)
    uploaded = st.file_uploader("Or upload a text gene list", type=["txt", "csv", "tsv"])
    if uploaded:
        pasted += "\n" + uploaded.getvalue().decode("utf-8")
    st.session_state.genes = sorted({
        token.strip().upper()
        for token in pasted.replace(",", " ").split()
        if token.strip()
    })
    st.write(f"{len(st.session_state.genes)} unique genes selected.")
else:
    st.info("All reference autosomes and comparison-valid sex chromosomes will be considered.")
