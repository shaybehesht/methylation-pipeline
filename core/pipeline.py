"""Synchronous pipeline orchestrator used by Streamlit."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import pandas as pd

from core.analysis import feasibility, intersect_and_rank, phenotype_segregation_rank
from core.annotations import genes_from_gtf, write_bed
from core.config import Affection, TrioConfig
from core.dmr import add_site_counts, read_regions, run_pair
from core.pileup import run_pileup
from core.plotting import effect_plot
from core.reasoning import explain, write_html_report
from core.thresholds import validate

Progress = Callable[[float, str], None]


def _scope_bed(config: TrioConfig, output: Path, gtf: str | None, cpg_islands: str | None) -> Path:
    bed = output / "scope.bed"
    if config.regions.mode == "targeted":
        if not gtf:
            raise ValueError("A GTF is required for targeted analysis")
        frame = genes_from_gtf(
            gtf, config.regions.genes,
            int(config.thresholds["promoter_pad"]), int(config.thresholds["body_pad"]),
        )
        missing = sorted(set(config.regions.genes) - set(frame["gene"]))
        if missing:
            raise ValueError(f"Genes absent from annotation: {', '.join(missing)}")
        return write_bed(frame, bed)
    if not cpg_islands:
        raise ValueError("A CpG-island BED is required for whole-genome and chromosome modes")
    islands = pd.read_csv(cpg_islands, sep="\t", names=["chrom", "start", "end", "name"])
    if config.regions.mode == "chromosomes":
        islands = islands[islands["chrom"].isin(config.regions.chromosomes)]
        if islands.empty:
            raise ValueError("No CpG islands found on the selected chromosomes")
    islands.to_csv(bed, sep="\t", header=False, index=False)
    return bed


def run(
    config: TrioConfig,
    gtf: str | None = None,
    cpg_islands: str | None = None,
    progress: Progress | None = None,
) -> dict:
    notify = progress or (lambda fraction, message: None)
    config.thresholds = validate(config.thresholds)
    output = config.ensure_output_dir()
    (output / "config.json").write_text(json.dumps(config.to_dict(), indent=2), encoding="utf-8")
    scope = _scope_bed(config, output, gtf, cpg_islands)
    log_path = output / "pipeline.log"
    pileups: dict[str, Path] = {}
    notify(0.05, "Preparing methylation pileups")
    with log_path.open("w", encoding="utf-8") as log:
        for index, sample in enumerate(config.samples):
            pileup_path = output / f"{sample.label}.bed.gz"
            pileups[sample.label] = run_pileup(
                sample.bam_path, str(pileup_path), config.reference_fasta,
                filter_threshold=float(config.thresholds["filter_threshold"]),
                include_bed=str(scope), log=log,
            )
            notify(0.1 + 0.15 * (index + 1), f"Pileup complete: {sample.label}")
        tables = []
        comparisons = config.comparisons()
        for index, comparison in enumerate(comparisons):
            path = output / f"{comparison.name}.tsv"
            run_pair(
                str(pileups[comparison.left.label]), str(pileups[comparison.right.label]),
                str(path), config.reference_fasta, log=log, regions=str(scope),
                min_valid_coverage=int(config.thresholds["min_valid_coverage"]),
            )
            table = add_site_counts(read_regions(path), pileups[comparison.left.label])
            table = table[table["chrom"].map(comparison.valid_chromosome)]
            tables.append(table)
            notify(0.6 + 0.08 * (index + 1), f"DMR complete: {comparison.name}")

    design = config.analysis_design()
    if design == "no_unaffected_control":
        raise ValueError(
            "Both relatives are marked affected. Add an unaffected comparator or "
            "set uncertain clinical status to unknown/not provided."
        )
    if design == "phenotype_segregation":
        r1, r2 = config.relatives
        if r1.affection == Affection.AFFECTED:
            similar, p_diff, affected_diff = tables[0], tables[1], tables[2]
            validators = (comparisons[0], comparisons[1], comparisons[2])
        else:
            similar, p_diff, affected_diff = tables[1], tables[0], tables[2].copy()
            affected_diff["effect"] = -affected_diff["effect"]
            validators = (comparisons[1], comparisons[0], comparisons[2])
        cutoff = float(config.thresholds["targeted_min_delta"]) / 100
        candidates, denominator_count = phenotype_segregation_rank(
            similar, p_diff, affected_diff,
            min_delta=cutoff,
            min_sites=int(config.thresholds["min_sites"]),
            max_pval=float(config.thresholds["max_pval"]),
            chromosome_validator=lambda chrom: all(
                comparison.valid_chromosome(chrom) for comparison in validators
            ),
        )
        denominator_label = "phenotype-discordant regions"
    else:
        candidates, cutoff = intersect_and_rank(
            tables[0], tables[1], tables[2],
            null_percentile=float(config.thresholds["null_percentile"]),
            min_sites=int(config.thresholds["min_sites"]),
            max_pval=float(config.thresholds["max_pval"]),
            chromosome_validator=lambda chrom: comparisons[0].valid_chromosome(chrom)
            and comparisons[1].valid_chromosome(chrom),
        )
        denominator_count = len(tables[2])
        denominator_label = "relative-null regions"
    if config.regions.mode == "targeted" and not candidates.empty:
        min_delta = float(config.thresholds["targeted_min_delta"]) / 100
        alpha = float(config.thresholds["alpha"])
        candidates = candidates[
            (candidates["mean_abs_effect"] >= min_delta)
            & (candidates["max_pvalue"] <= alpha)
        ].reset_index(drop=True)
        candidates["rank"] = range(1, len(candidates) + 1)
    candidates.to_csv(output / "proband_specific_DMRs.tsv", sep="\t", index=False)
    summary = feasibility(len(candidates), denominator_count)
    summary.update({"design": design, "denominator_label": denominator_label})
    reasoning = explain(summary, cutoff, config.caveats())
    figure = effect_plot(candidates, output / "dmr_effects.png")
    report = write_html_report(output / "report.html", "Methylation Trio Report", summary, reasoning, candidates, figure)
    result = {
        **summary, "null_cutoff": cutoff, "reasoning": reasoning,
        "evidence_status": config.evidence_status(),
        "report": str(report), "output": str(output),
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    notify(1.0, "Analysis complete")
    return result
