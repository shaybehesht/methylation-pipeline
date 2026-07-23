"""modkit pileup command construction, 5hmC folding, and bedMethyl indexing.

Mirrors the reference shell pipeline (scripts 01/01b/06/10):

* ``modkit pileup --cpg --combine-strands --modified-bases 5mC --filter-threshold``
* fold confident 5hmC calls (``N_other``, column 14) into ``N_canonical``
  (column 13) so ``valid_coverage == N_mod + N_canonical`` and ``modkit dmr
  pair`` accepts every row in single-site mode;
* sort, bgzip, and tabix-index the result.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pysam

from core.bedmethyl import bedmethyl_has_rows, pileup_modified_bases, validate_bedmethyl
from core.subprocess_util import run_checked


def build_command(
    bam: str,
    output: str,
    reference: str,
    *,
    filter_threshold: float = 0.7,
    region: str | None = None,
    combine_strands: bool = True,
    modified_bases: tuple[str, ...] = ("5mC",),
    threads: int | None = None,
    interval_size: int | None = None,
    suppress_progress: bool = True,
    log_filepath: str | None = None,
) -> list[str]:
    command = ["modkit", "pileup", bam, output, "--ref", reference]
    if region:
        command.extend(["--region", region])
    command.append("--cpg")
    if combine_strands:
        command.append("--combine-strands")
    # modkit >=0.6 requires an explicit modified base with --cpg; --filter-threshold
    # follows so the variadic modified-bases list terminates cleanly.
    command.append("--modified-bases")
    command.extend(modified_bases)
    command.extend(["--filter-threshold", str(filter_threshold)])
    if threads:
        command.extend(["--threads", str(threads)])
    if interval_size:
        command.extend(["--interval-size", str(interval_size)])
    if suppress_progress:
        command.append("--suppress-progress")
    if log_filepath:
        command.extend(["--log-filepath", log_filepath])
    return command


def fold_hmc(raw_line: str) -> str:
    """Fold column 14 (N_other, i.e. confident 5hmC) into column 13 (N_canonical).

    A confident 5hmC call is "not 5mC", so for a 5mC analysis it counts as
    unmodified. valid_coverage (col 10) and percent-modified (col 11) are
    unchanged; only the split between canonical and other moves.
    """

    fields = raw_line.rstrip("\n").split("\t")
    if len(fields) >= 14:
        try:
            fields[12] = str(int(fields[12]) + int(fields[13]))
            fields[13] = "0"
        except ValueError:
            return raw_line.rstrip("\n")
    return "\t".join(fields)


def _sort_bed(source: Path, dest: Path) -> None:
    """Coordinate-sort a BED file, preferring the system sorter for low memory."""

    if shutil.which("sort"):
        with dest.open("w") as handle:
            completed = subprocess.run(
                ["sort", "-k1,1", "-k2,2n", str(source)],
                stdout=handle, text=True,
                env={"LC_ALL": "C"},
            )
        if completed.returncode == 0:
            return
    frame = pd.read_csv(source, sep="\t", header=None, dtype=str)
    frame[1] = frame[1].astype(int)
    frame = frame.sort_values([0, 1]).astype({1: str})
    frame.to_csv(dest, sep="\t", header=False, index=False)


def fold_and_index(raw_bed: str | Path, output_gz: str | Path) -> Path:
    """Fold 5hmC, sort, bgzip, and tabix-index ``raw_bed`` into ``output_gz``."""

    raw_path = Path(raw_bed)
    out_path = Path(output_gz)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as scratch:
        folded = Path(scratch) / "folded.bed"
        with raw_path.open("r") as source, folded.open("w") as target:
            for line in source:
                if not line.strip():
                    continue
                target.write(fold_hmc(line) + "\n")
        ordered = Path(scratch) / "sorted.bed"
        _sort_bed(folded, ordered)
        if out_path.exists():
            out_path.unlink()
        pysam.tabix_compress(str(ordered), str(out_path), force=True)
    pysam.tabix_index(str(out_path), preset="bed", force=True)
    return out_path


def run_pileup(
    bam: str,
    output: str,
    reference: str,
    *,
    filter_threshold: float = 0.7,
    region: str | None = None,
    combine_strands: bool = True,
    modified_bases: tuple[str, ...] = ("5mC",),
    threads: int | None = None,
    include_bed: str | None = None,
    log=None,
) -> Path:
    """Pileup a BAM, fold 5hmC, and produce a bgzipped + tabixed bedMethyl.

    ``include_bed`` (a BED3 file) is honoured only as a modkit ``--include-bed``
    restriction; regions carrying names are handled by the targeted reader, not
    passed to modkit.
    """

    output_path = Path(output)
    raw_path = output_path.with_suffix(output_path.suffix + ".raw.bed")
    tabulated = pileup_modified_bases(modified_bases)
    modkit_threads = threads or 4
    with tempfile.TemporaryDirectory() as scratch:
        log_path = Path(scratch) / "modkit_pileup.log"
        command = build_command(
            bam, str(raw_path), reference,
            filter_threshold=filter_threshold, region=region,
            combine_strands=combine_strands, modified_bases=tabulated,
            threads=modkit_threads,
            interval_size=100_000,
            log_filepath=str(log_path),
        )
        if include_bed:
            command.extend(["--include-bed", include_bed])
        run_checked(command, log=log, log_filepath=log_path, tool="modkit pileup")
    try:
        fold_and_index(raw_path, output_path)
        n_rows, n_invalid, examples = validate_bedmethyl(output_path)
        if n_invalid:
            output_path.unlink(missing_ok=True)
            detail = "\n".join(examples[:5])
            raise RuntimeError(
                f"bedMethyl failed validation ({n_invalid}/{n_rows} bad rows) in {output_path}.\n"
                f"Examples:\n{detail}"
            )
    finally:
        raw_path.unlink(missing_ok=True)
    return output_path


def normalize_counts(frame: pd.DataFrame) -> pd.DataFrame:
    """Recompute N_other so valid, modified, canonical, and other sum consistently."""
    result = frame.copy()
    aliases = {
        "valid_coverage": ["valid_coverage", "N_valid_cov"],
        "modified": ["modified", "N_mod"],
        "canonical": ["canonical", "N_canonical"],
        "other": ["other", "N_other"],
    }
    resolved = {}
    for canonical, options in aliases.items():
        resolved[canonical] = next((name for name in options if name in result), None)
    if not all(resolved.values()):
        raise ValueError("Pileup must contain valid coverage, modified, canonical, and other counts")
    valid, modified, canonical, other = (resolved[k] for k in aliases)
    result[other] = (result[valid] - result[modified] - result[canonical]).clip(lower=0)
    return result
