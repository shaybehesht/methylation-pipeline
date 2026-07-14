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


def _run_segmentation(config: TrioConfig, output: Path, log, notify: Progress) -> dict:
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
    return {
        "candidates": candidates, "cutoff": cutoff,
        "denominator": denominator, "denominator_label": "relative-null segments",
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
    return {
        "candidates": candidates, "cutoff": float(config.thresholds["targeted_min_delta"]) / 100.0,
        "denominator": denominator, "denominator_label": "tested panel regions",
    }


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
        if config.regions.mode == "targeted":
            outcome = _run_targeted(config, output, gtf, log, notify)
        else:
            outcome = _run_segmentation(config, output, log, notify)

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
    report = write_html_report(
        output / "report.html", "Methylation Trio Report", summary, reasoning, candidates, figure
    )
    result = {
        **summary, "null_cutoff": cutoff, "reasoning": reasoning,
        "evidence_status": config.evidence_status(),
        "report": str(report), "output": str(output),
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    notify(1.0, "Analysis complete")
    return result
