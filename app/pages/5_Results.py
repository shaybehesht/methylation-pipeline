import json
from pathlib import Path

import pandas as pd
import streamlit as st

from app.state import initialize
from core.figures import ROLE_COLORS
from core.profile import parse_region, region_profile
from core.reporting import archive_size, build_run_archive, human_size

LARGE_ARCHIVE_BYTES = 500 * 1024 * 1024


def _run_pileups(run_dir: Path) -> tuple[dict, list]:
    """Return ({label: pileup_path}, ordered_labels) for a run, proband first."""
    labels: list[str] = []
    config_path = run_dir / "config.json"
    if config_path.exists():
        try:
            samples = json.loads(config_path.read_text(encoding="utf-8")).get("samples", [])
            proband = [s["label"] for s in samples if s.get("role") == "proband"]
            others = [s["label"] for s in samples if s.get("role") != "proband"]
            labels = proband + others
        except (ValueError, KeyError):
            labels = []
    if not labels:
        labels = sorted(p.name[:-7] for p in run_dir.glob("*.bed.gz"))
    pileups = {}
    for label in labels:
        path = run_dir / f"{label}.bed.gz"
        if path.exists() and (path.with_suffix(path.suffix + ".tbi")).exists():
            pileups[label] = path
    return pileups, [label for label in labels if label in pileups]

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

figdir = output / "figures"
gene_plots = sorted(figdir.glob("gene_*.png")) if figdir.exists() else []
overview_plots = []
if figdir.exists():
    for pattern in (
        "threshold_sweep.png", "wgs_karyotype.png", "effect_histogram.png",
        "targeted_heatmap_*.png",
    ):
        overview_plots.extend(sorted(figdir.glob(pattern)))

if gene_plots or overview_plots:
    st.subheader("Figures")
if gene_plots:
    names = [p.stem.replace("gene_", "") for p in gene_plots]
    chosen = st.selectbox("Per-gene methylation plot", names)
    st.image(str(figdir / f"gene_{chosen}.png"), use_container_width=True)
for plot in overview_plots:
    st.image(str(plot), caption=plot.stem.replace("_", " "), use_container_width=True)

figure = output / "dmr_effects.png"
if figure.exists():
    with st.expander("Summary effect scatter"):
        st.image(str(figure))

# ----------------------------------------------------------------- interactive
st.subheader("Interactive methylation profile (zoomable)")
pileups, ordered_labels = _run_pileups(output)
if not pileups:
    st.caption("Per-sample pileups were not found in this run directory.")
else:
    region_choices: dict[str, tuple[str, int, int]] = {}
    regions_tsv = output / "regions" / "dmr_regions.tsv"
    if regions_tsv.exists():
        try:
            reg = pd.read_csv(regions_tsv, sep="\t")
            for gene, group in reg[reg["region"] == "body"].groupby("gene"):
                row = group.iloc[0]
                region_choices[f"{gene} (gene)"] = (str(row["chrom"]), int(row["start"]), int(row["end"]))
        except Exception:
            pass
    dmr_tsv = output / "proband_specific_DMRs.tsv"
    if dmr_tsv.exists():
        try:
            dmrs = pd.read_csv(dmr_tsv, sep="\t")
            if {"chrom", "start", "end"} <= set(dmrs.columns):
                for _, row in dmrs.head(50).iterrows():
                    tag = f"DMR {row['chrom']}:{int(row['start'])}-{int(row['end'])}"
                    region_choices[tag] = (str(row["chrom"]), int(row["start"]), int(row["end"]))
        except Exception:
            pass

    options = list(region_choices) + ["Custom region…"]
    picked = st.selectbox("Region", options) if options else "Custom region…"
    if picked == "Custom region…":
        manual = st.text_input("Region (chrom:start-end)", placeholder="chr1:156,590,755-156,594,876")
        parsed = parse_region(manual) if manual else None
    else:
        parsed = region_choices[picked]

    controls = st.columns(2)
    flank = controls[0].slider("Flank (bp)", 0, 20000, 2000, 500)
    window = controls[1].slider("Sliding window (CpGs)", 1, 200, 20)

    if parsed:
        chrom, region_start, region_end = parsed
        view_start = max(0, region_start - flank)
        view_end = region_end + flank
        profile = region_profile(pileups, chrom, view_start, view_end, window=window)
        if profile.empty:
            st.info("No covered CpGs in this region for the selected samples.")
        else:
            try:
                import plotly.graph_objects as go

                colors = {label: ROLE_COLORS[i % len(ROLE_COLORS)] for i, label in enumerate(ordered_labels)}
                fig = go.Figure()
                for label in ordered_labels:
                    grp = profile[profile["sample"] == label]
                    if grp.empty:
                        continue
                    fig.add_trace(go.Scatter(
                        x=grp["pos"], y=grp["smooth"], mode="lines", name=label,
                        line=dict(color=colors.get(label), width=2),
                        hovertemplate="%{x:,} bp<br>%{y:.1f}%<extra>" + label + "</extra>",
                    ))
                fig.add_vrect(
                    x0=region_start, x1=region_end, fillcolor="#B0B0B0", opacity=0.18,
                    line_width=0, annotation_text="region", annotation_position="top left",
                )
                fig.update_layout(
                    height=520, hovermode="x unified",
                    xaxis_title="Genomic position", yaxis_title="mean % CpG methylation",
                    yaxis_range=[-3, 103], margin=dict(l=10, r=10, t=30, b=10),
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
                )
                fig.update_xaxes(rangeslider_visible=True)
                st.plotly_chart(fig, use_container_width=True)
                st.caption(
                    "Drag to zoom, double-click to reset, and use the mini-track below the "
                    "plot to pan across the region. The shaded band is the selected region."
                )
            except ModuleNotFoundError:
                st.warning("Install plotly for the interactive view: `pip install -e .` (adds plotly).")
    else:
        st.caption("Pick a region above or enter one as chrom:start-end.")
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
    "Bundle the entire run directory — manifests, logs, pileups, pairwise "
    "outputs, tables, figures, and the HTML report — into a single ZIP."
)
archive_path = output / "complete_run.zip"
if st.button("Build complete ZIP", type="primary"):
    with st.spinner("Packaging the complete run directory..."):
        archive_path = build_run_archive(output)
    st.success("Archive built.")

if archive_path.exists():
    size = archive_size(archive_path)
    st.write(f"Local archive: `{archive_path}` ({human_size(size)})")
    if size >= LARGE_ARCHIVE_BYTES:
        st.warning(
            f"This archive is large ({human_size(size)}). Downloading through the "
            "browser may be slow; the local path above can be copied directly."
        )
    st.download_button(
        "Download complete ZIP", archive_path.read_bytes(),
        file_name="complete_run.zip", mime="application/zip",
    )
