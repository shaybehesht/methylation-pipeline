"""Chromosome name resolution between references, BAMs, and bedMethyl."""
from __future__ import annotations

import gzip
from pathlib import Path

import pysam

from core.qc import fasta_contigs


def _chrom_variants(chrom: str) -> list[str]:
    stripped = chrom.removeprefix("chr")
    return list(dict.fromkeys([chrom, stripped, f"chr{stripped}"]))


def resolve_bam_contig(chrom: str, bam_paths: list[str]) -> str:
    """Return the contig name used in the modBAM header for ``chrom``."""
    for path in bam_paths:
        try:
            with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
                refs = set(bam.references)
        except (OSError, ValueError):
            continue
        for candidate in _chrom_variants(chrom):
            if candidate in refs:
                return candidate
    raise ValueError(
        f"Chromosome {chrom!r} was not found in the modBAM headers "
        f"(tried {_chrom_variants(chrom)})."
    )


def resolve_reference_contig(reference: str, chrom: str, bam_paths: list[str] | None = None) -> str:
    """Return the contig name present in ``reference`` for a logical chromosome."""
    ref_contigs = fasta_contigs(reference)
    for candidate in _chrom_variants(chrom):
        if candidate in ref_contigs:
            return candidate
    if bam_paths:
        for path in bam_paths:
            try:
                with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
                    bam_contigs = set(bam.references)
            except (OSError, ValueError):
                continue
            for candidate in _chrom_variants(chrom):
                if candidate in ref_contigs or candidate in bam_contigs:
                    for ref_name in _chrom_variants(candidate):
                        if ref_name in ref_contigs:
                            return ref_name
    raise ValueError(
        f"Chromosome {chrom!r} was not found in the reference "
        f"(tried {_chrom_variants(chrom)})."
    )


def pileup_contig_name(path: str | Path) -> str | None:
    """Return the chromosome/contig name used in a bedMethyl pileup file."""
    source = Path(path)
    if not source.exists() or source.stat().st_size == 0:
        return None
    opener = gzip.open if str(source).endswith(".gz") else open
    with opener(source, "rt") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            return line.split("\t", 1)[0]
    return None


def dmr_reference_contig(
    reference: str,
    chrom: str,
    pileup_paths: list[str],
    bam_paths: list[str] | None = None,
) -> tuple[str, str]:
    """Return ``(reference_contig, fasta_header)`` for modkit DMR on one chromosome.

    ``reference_contig`` selects sequence from the full FASTA. ``fasta_header`` is
    written into the single-chromosome FASTA so it matches bedMethyl row names.
    """
    ref_contig = resolve_reference_contig(reference, chrom, bam_paths)
    observed = [name for path in pileup_paths if (name := pileup_contig_name(path))]
    if observed:
        unique = set(observed)
        if len(unique) > 1:
            raise ValueError(
                f"Pileup contig names disagree for {chrom}: {sorted(unique)}"
            )
        return ref_contig, observed[0]
    return ref_contig, ref_contig
