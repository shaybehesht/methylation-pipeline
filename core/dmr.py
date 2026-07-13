"""modkit pairwise DMR wrapper and normalized table reader."""
from __future__ import annotations

import subprocess
import re
from pathlib import Path

import pandas as pd
import pysam
from scipy.stats import fisher_exact


def build_command(
    left: str,
    right: str,
    output: str,
    reference: str,
    *,
    regions: str | None = None,
    min_valid_coverage: int = 10,
) -> list[str]:
    command = [
        "modkit", "dmr", "pair", "-a", left, "-b", right, "-o", output,
        "--ref", reference, "--base", "C", "--header",
        "--min-valid-coverage", str(min_valid_coverage),
    ]
    if regions:
        command.extend(["--regions-bed", regions])
    return command


def run_pair(left: str, right: str, output: str, reference: str, log=None, **kwargs) -> Path:
    command = build_command(left, right, output, reference, **kwargs)
    completed = subprocess.run(command, text=True, stdout=log, stderr=log, check=False)
    if completed.returncode:
        raise RuntimeError(f"modkit dmr pair failed with exit code {completed.returncode}")
    return Path(output)


def read_regions(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [column.removeprefix("#").strip() for column in frame.columns]
    aliases = {
        "chrom": ("chrom", "chromosome", "contig"),
        "start": ("start", "region_start"),
        "end": ("end", "region_end"),
        "effect": ("effect", "effect_size", "delta"),
        "pvalue": ("pvalue", "p_value", "map_p_value"),
        "n_sites": ("n_sites", "N_sites", "num_sites", "site_count"),
    }
    rename = {}
    for canonical, candidates in aliases.items():
        source = next((column for column in candidates if column in frame), None)
        if source is None and canonical not in {"pvalue", "n_sites"}:
            raise ValueError(f"Cannot find {canonical} in {path}")
        if source:
            rename[source] = canonical
    frame = frame.rename(columns=rename)
    if "n_sites" not in frame:
        name_column = next((name for name in ("name", "region_name") if name in frame), None)
        frame["n_sites"] = (
            frame[name_column].astype(str).map(
                lambda value: int(re.search(r"(?:CpG:\s*|n=)(\d+)", value).group(1))
                if re.search(r"(?:CpG:\s*|n=)(\d+)", value) else 1
            ) if name_column else 1
        )
    if "pvalue" not in frame:
        required = {"a_total", "b_total", "a_pct_modified", "b_pct_modified"}
        if not required <= set(frame):
            frame["pvalue"] = 1.0
        else:
            def region_pvalue(row) -> float:
                a_mod = round(float(row["a_total"]) * float(row["a_pct_modified"]))
                b_mod = round(float(row["b_total"]) * float(row["b_pct_modified"]))
                return float(fisher_exact([
                    [a_mod, int(row["a_total"]) - a_mod],
                    [b_mod, int(row["b_total"]) - b_mod],
                ]).pvalue)
            frame["pvalue"] = frame.apply(region_pvalue, axis=1)
    return frame


def add_site_counts(frame: pd.DataFrame, pileup_path: str | Path) -> pd.DataFrame:
    """Count covered CpG records per scored region from the indexed pileup."""
    result = frame.copy()
    with pysam.TabixFile(str(pileup_path)) as pileup:
        counts = []
        for row in result.itertuples():
            try:
                positions = {
                    int(line.split("\t", 3)[1])
                    for line in pileup.fetch(str(row.chrom), int(row.start), int(row.end))
                }
                counts.append(len(positions))
            except (ValueError, OSError):
                counts.append(0)
    result["n_sites"] = counts
    return result
