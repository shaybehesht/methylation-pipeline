"""BAM/reference compatibility and modified-read QC gates."""
from __future__ import annotations

from pathlib import Path

import pysam


def fasta_contigs(path: str) -> dict[str, int]:
    with pysam.FastaFile(path) as fasta:
        return dict(zip(fasta.references, fasta.lengths))


def has_mn_tag(path: str, reads_to_check: int = 1000) -> bool | None:
    """Whether the modBAM carries the ``MN`` tag needed for --combine-strands.

    Returns ``True``/``False`` when it can be determined from the first mapped
    reads, or ``None`` when the file cannot be read (so callers can avoid
    overriding user intent on an ambiguous result). ONT modBAMs (Dorado) carry
    ``MN``; PacBio HiFi modBAMs typically do not.
    """
    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
            checked = 0
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                checked += 1
                if read.has_tag("MN"):
                    return True
                if checked >= reads_to_check:
                    break
    except (OSError, ValueError):
        return None
    return False


def any_mn_tag(paths: list[str], reads_to_check: int = 1000) -> bool:
    """True if any BAM has an ``MN`` tag, or if none could be inspected.

    Used to decide whether to keep ``--combine-strands``: only when every BAM is
    readable *and* none carry ``MN`` (e.g. a PacBio trio) do we return ``False``.
    """
    determinable = False
    for path in paths:
        result = has_mn_tag(path, reads_to_check)
        if result is True:
            return True
        if result is False:
            determinable = True
    return not determinable


def inspect_bam(path: str, reference: str | None = None, reads_to_check: int = 2000) -> dict:
    if not Path(path).exists():
        raise FileNotFoundError(path)
    with pysam.AlignmentFile(path, "rb", check_sq=True) as bam:
        bam_contigs = dict(zip(bam.references, bam.lengths))
        header = bam.header.to_dict()
        hp_reads = checked = mn_reads = 0
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            checked += 1
            hp_reads += int(read.has_tag("HP"))
            mn_reads += int(read.has_tag("MN"))
            if checked >= reads_to_check:
                break
    reference_matches = None
    mismatch = []
    if reference:
        ref_contigs = fasta_contigs(reference)
        mismatch = [
            name for name, length in bam_contigs.items()
            if name not in ref_contigs or ref_contigs[name] != length
        ]
        reference_matches = not mismatch
    program_text = " ".join(
        f"{entry.get('PN', '')} {entry.get('VN', '')} {entry.get('CL', '')}"
        for entry in header.get("PG", [])
    )
    model = next(
        (token.split("=", 1)[1] for token in program_text.split() if "model=" in token.lower()),
        "not recorded in BAM header",
    )
    return {
        "bam": path,
        "reference_matches": reference_matches,
        "mismatched_contigs": mismatch[:20],
        "reads_checked": checked,
        "hp_fraction": hp_reads / checked if checked else 0.0,
        "has_hp_tags": hp_reads > 0,
        "has_mn_tags": mn_reads > 0,
        "basecaller_model": model,
    }


def gate(results: list[dict]) -> tuple[bool, list[str]]:
    errors = []
    models = {result["basecaller_model"] for result in results}
    for result in results:
        if result["reference_matches"] is False:
            errors.append(f"{result['bam']}: BAM contigs do not match the reference")
    if len(models) > 1:
        errors.append("BAMs report different basecaller models (possible batch effect)")
    return not errors, errors
