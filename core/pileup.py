"""modkit pileup command construction and bedMethyl normalization."""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd


def build_command(
    bam: str,
    output: str,
    reference: str,
    *,
    filter_threshold: float = 0.7,
    region: str | None = None,
    include_bed: str | None = None,
    modified_bases: tuple[str, ...] = ("5mC",),
) -> list[str]:
    command = [
        "modkit", "pileup", bam, output, "--ref", reference,
        "--filter-threshold", str(filter_threshold),
    ]
    # modkit >=0.6 requires --modified-bases whenever --cpg is used; it also
    # enables the optimized CpG routines. --cpg follows so the variadic list of
    # modified bases terminates cleanly.
    if modified_bases:
        command.append("--modified-bases")
        command.extend(modified_bases)
    command.extend(["--cpg", "--bgzf"])
    if region:
        command.extend(["--region", region])
    if include_bed:
        command.extend(["--include-bed", include_bed])
    return command


def run_pileup(*args, log=None, **kwargs) -> Path:
    command = build_command(*args, **kwargs)
    completed = subprocess.run(command, text=True, stdout=log, stderr=log, check=False)
    if completed.returncode:
        raise RuntimeError(f"modkit pileup failed with exit code {completed.returncode}")
    indexed = subprocess.run(["tabix", "-f", "-p", "bed", args[1]], text=True, stdout=log, stderr=log, check=False)
    if indexed.returncode:
        raise RuntimeError(f"tabix indexing failed with exit code {indexed.returncode}")
    return Path(args[1])


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
