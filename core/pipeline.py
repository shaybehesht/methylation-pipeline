"""Synchronous pipeline orchestrator used by Streamlit.

Two analysis modes mirror the reference scripts:

* whole-genome / selected chromosomes -> per-chromosome ``modkit dmr pair
  --segment`` (HMM segmentation, memory-bounded by a one-chromosome FASTA)
  followed by the empirical-null intersection of scripts 02/11;
* targeted -> per-promoter/gene-body paired Wilcoxon from the pileups, the
  standalone method of scripts 06/07.

Pileups always fold confident 5hmC into the canonical count so ``modkit dmr``
accepts every row (script 01b).
"""
from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from pathlib import Path
from typing import Callable

import pandas as pd
import pysam

from core.analysis import feasibility
from core.annotations import panel_regions, write_bed3
from core.config import Sex, TrioConfig
from core.dmr import read_segments, run_pair
from core.pileup import run_pileup
from core.plotting import effect_plot
from core.reasoning import explain, write_html_report
from core.segments import proband_private_dmrs, significant
from core.targeted import score_regions
from core.thresholds import validate

Progress = Callable[[float, str], None]

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
WHOLE_GENOME_CHROMS = AUTOSOMES + ["chrX"]


def _resolve_annotations(config: TrioConfig, gtf: str | None, cpg_islands: str | None):
    if (gtf is None or cpg_islands is None) and config.assembly:
        from core.references import prepared_paths

        paths = prepared_paths(config.assembly)
        gtf = gtf or str(paths["gtf"])
        cpg_islands = cpg_islands or str(paths["cpg_islands"])
    return gtf, cpg_islands


def _one_chromosome_fasta(reference: str, chrom: str, dest: Path) -> Path:
    """Write and index a single-chromosome FASTA to bound modkit's memory use."""
    if not dest.exists():
        if shutil.which("samtools"):
            with dest.open("w") as handle:
                subprocess.run(["samtools", "faidx", reference, chrom], stdout=handle, check=True)
        else:
            dest.write_text(pysam.faidx(reference, chrom))
        pysam.faidx(str(dest))
    return dest


def _merge_pileups(parts: list[Path], destination: Path) -> Path:
    """Concatenate coordinate-sorted per-chromosome pileups into one indexed file."""
    scratch = destination.with_suffix(destination.suffix + ".merge.bed")
    with scratch.open("w") as target:
        for part in parts:
            if not part.exists():
                continue
            with gzip.open(part, "rt") as source:
                shutil.copyfileobj(source, target)
    if destination.exists():
        destination.unlink()
    pysam.tabix_compress(str(scratch), str(destination), force=True)
    scratch.unlink(missing_ok=True)
    pysam.tabix_index(str(destination), preset="bed", force=True)
    return destination


def _run_segmentation(config: TrioConfig, output: Path, gtf: str | None, log, notify: Progress) -> dict:
    reference = config.reference_fasta
    comparisons = config.comparisons()
    samples = config.samples
    if config.regions.mode == "whole_genome":
        chroms = list(WHOLE_GENOME_CHROMS)
    else:
        chroms = list(config.regions.chromosomes)
    chroms = [c for c in chroms if c.lower() not in {"chrm", "chrmt", "chry"}]

    split = output / "split"
    dmr_dir = output / "dmr"
    for directory in (split, dmr_dir):
        directory.mkdir(parents=True, exist_ok=True)

    combine = config.combine_strands
    modified = tuple(config.modified_bases) or ("5mC",)
    min_cov = int(config.thresholds["min_valid_coverage"])
    total_units = max(len(chroms), 1)

    per_sample_parts: dict[str, list[Path]] = {sample.label: [] for sample in samples}
    for index, chrom in enumerate(chroms):
        notify(0.05 + 0.5 * index / total_units, f"Pileup + DMR: {chrom}")
        chrom_pileups: dict[str, Path] = {}
        for sample in samples:
            part = split / f"{sample.label}.{chrom}.bed.gz"
            if not part.exists():
                run_pileup(
                    sample.bam_path, str(part), reference,
                    filter_threshold=float(config.thresholds["filter_threshold"]),
                    region=chrom, combine_strands=combine, modified_bases=modified, log=log,
                )
            chrom_pileups[sample.label] = part
            per_sample_parts[sample.label].append(part)

        chrom_fasta = _one_chromosome_fasta(reference, chrom, split / f"{chrom}.fa")
        for comparison in comparisons:
            if not comparison.valid_chromosome(chrom):
                continue
            segment_path = dmr_dir / f"{comparison.name}.{chrom}.segments.bed"
            run_pair(
                str(chrom_pileups[comparison.left.label]),
                str(chrom_pileups[comparison.right.label]),
                str(dmr_dir / f"{comparison.name}.{chrom}.sites.bed"),
                str(chrom_fasta), log=log, segment=str(segment_path),
                min_valid_coverage=min_cov,
            )

    for stale in split.glob("*.fa"):
        stale.unlink(missing_ok=True)
    for stale in split.glob("*.fa.fai"):
        stale.unlink(missing_ok=True)

    notify(0.6, "Merging segments and pileups")
    empty_segments = pd.DataFrame(
        columns=["chrom", "start", "end", "state", "num_sites", "effect", "ci_excludes_zero"]
    )
    segment_frames: dict[str, pd.DataFrame] = {}
    for comparison in comparisons:
        frames = [
            read_segments(dmr_dir / f"{comparison.name}.{chrom}.segments.bed")
            for chrom in chroms
            if (dmr_dir / f"{comparison.name}.{chrom}.segments.bed").exists()
        ]
        merged = pd.concat(frames, ignore_index=True) if frames else empty_segments.copy()
        segment_frames[comparison.name] = merged
        merged.to_csv(output / f"{comparison.name}.segments.tsv", sep="\t", index=False)

    for sample in samples:
        _merge_pileups(per_sample_parts[sample.label], output / f"{sample.label}.bed.gz")

    pm = segment_frames[comparisons[0].name]
    pb = segment_frames[comparisons[1].name]
    mb = segment_frames[comparisons[2].name]

    def validator(chrom: str) -> bool:
        return comparisons[0].valid_chromosome(chrom) and comparisons[1].valid_chromosome(chrom)

    candidates, cutoff = proband_private_dmrs(
        pm, pb, mb,
        null_percentile=float(config.thresholds["null_percentile"]),
        min_sites=int(config.thresholds["min_sites"]),
        require_ci=True,
        chromosome_validator=validator,
    )
    if not candidates.empty:
        candidates["mean_abs_effect"] = candidates[["effect_1", "effect_2"]].abs().mean(axis=1)
    null_variable = significant(
        mb, cutoff, int(config.thresholds["min_sites"]), require_ci=True
    )
    denominator = max(len(null_variable), 1)

    notify(0.95, "Rendering figures")
    from core import figures
    from core.annotations import all_genes, annotate_with_genes
    from core.config import Relationship
    from core.segments import flip_segments, private_count

    figdir = output / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # Annotate candidate DMRs with overlapping genes for the karyotype labels.
    if gtf and not candidates.empty:
        try:
            candidates = annotate_with_genes(candidates, all_genes(gtf))
        except (OSError, ValueError):
            pass

    # Threshold sweep: each individual's private-DMR count vs |effect| (autosomes).
    autosomes = set(AUTOSOMES)

    def autosomal(frame: pd.DataFrame) -> pd.DataFrame:
        return frame[frame["chrom"].isin(autosomes)] if not frame.empty else frame

    proband, (relative_one, relative_two) = config.proband, config.relatives
    pm_a, pb_a, mb_a = autosomal(pm), autosomal(pb), autosomal(mb)
    focal = {
        proband.label: (pm_a, pb_a, mb_a),
        relative_one.label: (flip_segments(pm_a), mb_a, pb_a),
        relative_two.label: (flip_segments(pb_a), flip_segments(mb_a), pm_a),
    }
    thresholds = [round(0.10 + 0.025 * step, 4) for step in range(24)]
    min_sites = int(config.thresholds["min_sites"])
    counts = {
        label: [private_count(c1, c2, null, t, min_sites) for t in thresholds]
        for label, (c1, c2, null) in focal.items()
    }

    def role_style(sample) -> tuple[str, str, str]:
        if sample.role.value == "proband":
            return figures.ROLE_COLORS[0], "-", f"{sample.label} (proband)"
        if sample.relationship == Relationship.SIBLING:
            return figures.ROLE_COLORS[2], "-", f"{sample.label} (sibling control)"
        if sample.relationship in (Relationship.MOTHER, Relationship.FATHER):
            return figures.ROLE_COLORS[1], "--", f"{sample.label} (parent)"
        return "#7F8C8D", "--", sample.label

    series = []
    for sample in (proband, relative_one, relative_two):
        color, linestyle, label = role_style(sample)
        series.append({
            "label": label, "counts": counts[sample.label], "color": color,
            "linestyle": linestyle, "linewidth": 2.4 if linestyle == "-" else 1.4,
            "alpha": 1.0 if linestyle == "-" else 0.6,
        })

    produced = [
        figures.sweep_plot(
            thresholds, series, figdir / "threshold_sweep.png",
            ylabel="private DMRs (autosomes)",
            title="Proband vs the family across effect thresholds",
        ),
        figures.karyotype_plot(
            candidates, figdir / "wgs_karyotype.png",
            title=f"Proband-private DMRs, |effect| >= {cutoff:.3g}",
        ),
        figures.effect_histogram(candidates, mb, figdir / "effect_histogram.png"),
    ]
    return {
        "candidates": candidates, "cutoff": cutoff,
        "denominator": denominator, "denominator_label": "relative-null segments",
        "figures": produced,
    }


def _run_targeted(config: TrioConfig, output: Path, gtf: str | None, log, notify: Progress) -> dict:
    if not gtf:
        raise ValueError("A GENCODE GTF is required for targeted analysis")
    reference = config.reference_fasta
    named, extract, missing = panel_regions(
        gtf, config.regions.genes,
        int(config.thresholds["promoter_pad"]), int(config.thresholds["body_pad"]),
    )
    if named.empty:
        raise ValueError(f"No panel genes found in the annotation: {', '.join(config.regions.genes)}")
    if missing:
        log.write(f"Genes not found in annotation: {', '.join(missing)}\n")
    regions_dir = output / "regions"
    regions_dir.mkdir(parents=True, exist_ok=True)
    named.to_csv(regions_dir / "dmr_regions.tsv", sep="\t", index=False)
    extract_bed = write_bed3(extract, regions_dir / "extract.bed")

    combine = config.combine_strands
    modified = tuple(config.modified_bases) or ("5mC",)
    pileups: dict[str, Path] = {}
    for index, sample in enumerate(config.samples):
        notify(0.1 + 0.25 * (index + 1), f"Pileup (panel): {sample.label}")
        path = output / f"{sample.label}.bed.gz"
        run_pileup(
            sample.bam_path, str(path), reference,
            filter_threshold=float(config.thresholds["filter_threshold"]),
            combine_strands=combine, modified_bases=modified,
            include_bed=str(extract_bed), log=log,
        )
        pileups[sample.label] = path

    proband = config.proband
    relative_one, relative_two = config.relatives
    female = next((r.label for r in config.relatives if r.sex == Sex.FEMALE), None)
    notify(0.8, "Scoring panel regions")
    scored = score_regions(
        named, pileups,
        proband=proband.label, relative_one=relative_one.label, relative_two=relative_two.label,
        female_relative=female,
        min_cov=int(config.thresholds["min_valid_coverage"]),
        min_cpgs=int(config.thresholds["min_sites"]),
        min_delta=float(config.thresholds["targeted_min_delta"]),
        alpha=float(config.thresholds["alpha"]),
    )
    scored.to_csv(output / "targeted_results.tsv", sep="\t", index=False)
    candidates = scored[scored.get("candidate", False)].copy() if not scored.empty else scored
    if not candidates.empty:
        candidates = candidates.reset_index(drop=True)
        candidates.insert(0, "rank", range(1, len(candidates) + 1))
        candidates["effect_1"] = candidates["delta_p_r1"] / 100.0
        candidates["mean_abs_effect"] = candidates["min_abs_delta"] / 100.0
    denominator = max(int((scored["n_cpgs"] >= int(config.thresholds["min_sites"])).sum()), 1) if not scored.empty else 1

    notify(0.9, "Rendering figures")
    figure_paths = _targeted_figures(
        output, named, scored, pileups, proband, relative_one, relative_two,
        int(config.thresholds["min_valid_coverage"]),
    )
    return {
        "candidates": candidates, "cutoff": float(config.thresholds["targeted_min_delta"]) / 100.0,
        "denominator": denominator, "denominator_label": "tested panel regions",
        "figures": figure_paths,
    }


def _targeted_figures(output, named, scored, pileups, proband, relative_one, relative_two, min_cov):
    """Per-gene methylation plots and promoter/body heatmaps (mirrors script 07)."""
    from core import figures

    figdir = output / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    labels = [proband.label, relative_one.label, relative_two.label]
    colors = {label: figures.ROLE_COLORS[index] for index, label in enumerate(labels)}
    pileup_order = {label: pileups[label] for label in labels}
    meth_columns = ["proband_meth", f"{relative_one.label}_meth", f"{relative_two.label}_meth"]

    produced: list = []
    if not scored.empty and all(column in scored.columns for column in meth_columns):
        for kind in ("promoter", "body"):
            path = figures.targeted_heatmap(
                scored, kind, meth_columns, labels, figdir / f"targeted_heatmap_{kind}.png"
            )
            if path:
                produced.append(path)

    for gene, group in list(named.groupby("gene"))[:80]:
        body = group[group["region"] == "body"]
        promoter = group[group["region"] == "promoter"]
        if body.empty:
            continue
        body_row = body.iloc[0]
        promoter_row = promoter.iloc[0] if not promoter.empty else body_row
        subtitle = ""
        if not scored.empty:
            match = scored[(scored["gene"] == gene) & (scored["region"] == "promoter")]
            if not match.empty and pd.notna(match.iloc[0].get(meth_columns[0])):
                row = match.iloc[0]
                subtitle = "promoter " + "  ".join(
                    f"{label[0].upper()} {row[column]:.0f}%"
                    for label, column in zip(labels, meth_columns)
                )
        path = figures.gene_locus_plot(
            gene, str(body_row["chrom"]), int(body_row["start"]), int(body_row["end"]),
            (int(promoter_row["start"]), int(promoter_row["end"])),
            pileup_order, figdir / f"gene_{gene}.png",
            min_cov=min_cov, subtitle=subtitle, colors=colors,
        )
        if path:
            produced.append(path)
    return produced


def run(
    config: TrioConfig,
    gtf: str | None = None,
    cpg_islands: str | None = None,
    progress: Progress | None = None,
) -> dict:
    notify = progress or (lambda fraction, message: None)
    gtf, cpg_islands = _resolve_annotations(config, gtf, cpg_islands)
    config.thresholds = validate(config.thresholds)
    output = config.ensure_output_dir()
    (output / "config.json").write_text(json.dumps(config.to_dict(), indent=2), encoding="utf-8")
    log_path = output / "pipeline.log"
    notify(0.02, "Starting analysis")
    with log_path.open("w", encoding="utf-8") as log:
        from core.bam_index import ensure_index

        for sample in config.samples:
            ensure_index(sample.bam_path, log)
        # PacBio HiFi (and other) modBAMs without an MN tag cannot use modkit
        # --combine-strands; auto-disable it so those inputs "just work".
        if config.combine_strands:
            from core.qc import any_mn_tag

            if not any_mn_tag([sample.bam_path for sample in config.samples]):
                config.combine_strands = False
                log.write(
                    "No MN tags detected in the modBAMs (e.g. PacBio HiFi); "
                    "disabling --combine-strands automatically.\n"
                )
        if config.regions.mode == "targeted":
            outcome = _run_targeted(config, output, gtf, log, notify)
        else:
            outcome = _run_segmentation(config, output, gtf, log, notify)

    candidates = outcome["candidates"]
    cutoff = outcome["cutoff"]
    candidates.to_csv(output / "proband_specific_DMRs.tsv", sep="\t", index=False)
    if not candidates.empty and {"chrom", "start", "end"} <= set(candidates.columns):
        bed = candidates[["chrom", "start", "end"]].copy()
        bed["name"] = [f"DMR{i + 1}" for i in range(len(candidates))]
        bed.to_csv(output / "proband_specific_DMRs.bed", sep="\t", header=False, index=False)

    summary = feasibility(len(candidates), outcome["denominator"])
    summary.update({
        "design": config.analysis_design(),
        "denominator_label": outcome["denominator_label"],
    })
    reasoning = explain(summary, cutoff, config.caveats())
    figure = effect_plot(candidates, output / "dmr_effects.png")
    all_figures = list(outcome.get("figures", [])) + [figure]
    all_figures = [fig for fig in all_figures if fig is not None]
    report = write_html_report(
        output / "report.html", "MANGO — Methylation Analysis for Novel Genomic Outcomes",
        summary, reasoning, candidates, figures=all_figures,
    )
    result = {
        **summary, "null_cutoff": cutoff, "reasoning": reasoning,
        "evidence_status": config.evidence_status(),
        "report": str(report), "output": str(output),
        "figures": [str(fig) for fig in all_figures],
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    notify(1.0, "Analysis complete")
    return result
