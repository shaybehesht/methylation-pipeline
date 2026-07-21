"""Headless command-line entry point for the MANGO trio DMR pipeline.

The Streamlit app is one caller of :func:`core.pipeline.run`; this module is a
second, non-interactive caller so the exact same analysis can run in a batch
context such as a WDL task on AnVIL / Terra (Cromwell). It builds a validated
:class:`core.config.TrioConfig` from command-line flags and runs the pipeline,
writing every output under ``--output-dir``.

References are passed as explicit files (``--reference-fasta``, ``--gtf``,
``--cpg-islands``) so a cloud task is hermetic and needs no network; if the GTF
and CpG-island tracks are omitted, the managed-assembly download path in
``core.pipeline`` is used as a fallback.

Example::

    mango-run \\
        --proband-bam proband.bam --proband-sex F --proband-affection affected \\
        --relative1-bam mother.bam --relative1-sex F --relative1-relationship mother \\
        --relative2-bam father.bam --relative2-sex M --relative2-relationship father \\
        --reference-fasta hg38.fa --gtf gencode.gtf --cpg-islands cpg.txt \\
        --mode targeted --genes MECP2 UBE3A SNRPN \\
        --output-dir out
"""
from __future__ import annotations

import argparse
import sys
from typing import Sequence

from core.config import (
    Affection,
    RegionConfig,
    Relationship,
    Role,
    Sample,
    Sex,
    TrioConfig,
)
from core.thresholds import REGISTRY, defaults


def _add_sample_flags(parser: argparse.ArgumentParser, slot: str, *, proband: bool) -> None:
    group = parser.add_argument_group(f"{slot} sample")
    group.add_argument(f"--{slot}-bam", required=True, help=f"modBAM path for the {slot}.")
    group.add_argument(f"--{slot}-label", default=None, help=f"Display label (default: {slot}).")
    group.add_argument(
        f"--{slot}-sex", required=True, choices=[s.value for s in Sex],
        help="Sample sex (F/M); used to decide chrX/chrY comparison validity.",
    )
    group.add_argument(
        f"--{slot}-affection", default=None, choices=[a.value for a in Affection],
        help="Optional clinical status.",
    )
    if not proband:
        group.add_argument(
            f"--{slot}-relationship", default=None,
            choices=[r.value for r in Relationship],
            help="Optional relationship to the proband.",
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mango-run",
        description="Run the MANGO three-sample methylation DMR pipeline headlessly.",
    )
    _add_sample_flags(parser, "proband", proband=True)
    _add_sample_flags(parser, "relative1", proband=False)
    _add_sample_flags(parser, "relative2", proband=False)

    refs = parser.add_argument_group("reference")
    refs.add_argument(
        "--reference-fasta", default="",
        help="Reference FASTA (indexed). Omit to auto-download via --assembly.",
    )
    refs.add_argument("--gtf", default=None, help="GENCODE GTF (required for targeted mode).")
    refs.add_argument("--cpg-islands", default=None, help="UCSC CpG-island track.")
    refs.add_argument(
        "--assembly", default="",
        help="Managed assembly key (hg38/hg19). When --reference-fasta is omitted, the "
        "matching FASTA, GENCODE GTF, and CpG islands are downloaded and prepared "
        "automatically (cached in METHYL_TRIO_REFERENCE_CACHE).",
    )

    regions = parser.add_argument_group("regions")
    regions.add_argument(
        "--mode", default="targeted",
        choices=["whole_genome", "chromosomes", "targeted"],
        help="Analysis scope (default: targeted).",
    )
    regions.add_argument(
        "--chromosomes", nargs="+", default=None,
        help="Chromosomes for --mode chromosomes (e.g. chr1 chr2 chr11 chr15).",
    )
    regions.add_argument(
        "--genes", nargs="+", default=None,
        help="Gene symbols for --mode targeted (e.g. MECP2 UBE3A SNRPN).",
    )

    parser.add_argument("--phased-vcf", default="", help="Optional phased family VCF.")
    parser.add_argument(
        "--modified-base", action="append", dest="modified_bases", default=None,
        choices=["5mC", "5hmC"], help="Modified base(s) to tabulate (default: 5mC).",
    )
    parser.add_argument(
        "--no-combine-strands", action="store_true",
        help="Disable modkit --combine-strands (for modBAMs lacking MN tags).",
    )
    parser.add_argument(
        "--set-threshold", action="append", dest="thresholds", default=None,
        metavar="KEY=VALUE",
        help="Override a threshold; repeatable. Keys: " + ", ".join(REGISTRY),
    )
    parser.add_argument("--output-dir", required=True, help="Directory for all run outputs.")
    return parser


def parse_thresholds(pairs: Sequence[str] | None) -> dict[str, float | int]:
    """Parse ``KEY=VALUE`` override strings into a threshold dictionary."""
    result: dict[str, float | int] = defaults()
    for raw in pairs or []:
        if "=" not in raw:
            raise ValueError(f"Threshold override must be KEY=VALUE, got {raw!r}")
        key, _, value = raw.partition("=")
        key = key.strip()
        if key not in REGISTRY:
            raise ValueError(f"Unknown threshold {key!r}; choose from {', '.join(REGISTRY)}")
        text = value.strip()
        number: float | int = int(text) if text.lstrip("-").isdigit() else float(text)
        result[key] = number
    return result


def resolve_reference(args: argparse.Namespace, progress=None) -> None:
    """Ensure ``args.reference_fasta`` is set, auto-downloading from ``--assembly``.

    When no ``--reference-fasta`` was given, the managed assembly (hg38/hg19) is
    downloaded and prepared, and the FASTA, GTF, and CpG-island paths are filled
    in (any explicitly provided ``--gtf``/``--cpg-islands`` take precedence). This
    performs network I/O, so it is kept out of :func:`build_config` (which stays
    pure and unit-testable).
    """
    if args.reference_fasta:
        return
    if not args.assembly:
        raise ValueError(
            "Provide --reference-fasta, or --assembly (hg38/hg19) to auto-download a "
            "managed reference."
        )
    from core.references import prepare_assembly

    paths = prepare_assembly(args.assembly, progress=progress)
    args.reference_fasta = str(paths["fasta"])
    args.gtf = args.gtf or str(paths["gtf"])
    args.cpg_islands = args.cpg_islands or str(paths["cpg_islands"])


def _sample(args: argparse.Namespace, slot: str, role: Role) -> Sample:
    bam = getattr(args, f"{slot}_bam")
    label = getattr(args, f"{slot}_label") or slot
    sex = Sex(getattr(args, f"{slot}_sex"))
    affection_value = getattr(args, f"{slot}_affection")
    relationship_value = getattr(args, f"{slot}_relationship", None)
    return Sample(
        label=label,
        bam_path=bam,
        sex=sex,
        role=role,
        relationship=Relationship(relationship_value) if relationship_value else None,
        affection=Affection(affection_value) if affection_value else None,
    )


def build_config(args: argparse.Namespace) -> TrioConfig:
    """Build a validated :class:`TrioConfig` from parsed arguments."""
    samples = [
        _sample(args, "proband", Role.PROBAND),
        _sample(args, "relative1", Role.RELATIVE),
        _sample(args, "relative2", Role.RELATIVE),
    ]
    region_kwargs: dict[str, object] = {"mode": args.mode}
    if args.chromosomes:
        region_kwargs["chromosomes"] = list(args.chromosomes)
    if args.genes:
        region_kwargs["genes"] = list(args.genes)
    regions = RegionConfig(**region_kwargs)

    return TrioConfig(
        samples=samples,
        reference_fasta=args.reference_fasta,
        output_dir=args.output_dir,
        regions=regions,
        thresholds=parse_thresholds(args.thresholds),
        phased_vcf=args.phased_vcf or "",
        assembly=args.assembly or "",
        modified_bases=list(args.modified_bases) if args.modified_bases else ["5mC"],
        combine_strands=not args.no_combine_strands,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    def report(fraction: float, message: str) -> None:
        print(f"[{fraction * 100:5.1f}%] {message}", file=sys.stderr, flush=True)

    try:
        resolve_reference(args, progress=report)
        config = build_config(args)
    except (ValueError, KeyError) as exc:
        parser.error(str(exc))

    from core.pipeline import run

    result = run(config, gtf=args.gtf, cpg_islands=args.cpg_islands, progress=report)
    print(result.get("output", args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
