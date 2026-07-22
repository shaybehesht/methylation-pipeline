"""modkit pairwise DMR wrapper and normalized table reader."""
from __future__ import annotations

import subprocess
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact


def build_command(
    left: str,
    right: str,
    output: str,
    reference: str,
    *,
    segment: str | None = None,
    regions: str | None = None,
    min_valid_coverage: int = 10,
    missing: str | None = None,
    force: bool = True,
    threads: int | None = None,
    log_filepath: str | None = None,
) -> list[str]:
    command = [
        "modkit", "dmr", "pair", "-a", left, "-b", right, "-o", output,
        "--ref", reference, "--base", "C", "--header",
        "--min-valid-coverage", str(min_valid_coverage),
    ]
    if force:
        command.append("--force")
    if segment:
        command.extend(["--segment", segment])
    if regions:
        command.extend(["--regions-bed", regions])
    if missing:
        command.extend(["--missing", missing])
    if threads:
        command.extend(["--threads", str(threads)])
    if log_filepath:
        command.extend(["--log-filepath", log_filepath])
    return command


def run_pair(left: str, right: str, output: str, reference: str, log=None, **kwargs) -> Path:
    command = build_command(left, right, output, reference, **kwargs)
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if log is not None:
        if completed.stdout:
            log.write(completed.stdout)
        if completed.stderr:
            log.write(completed.stderr)
    if completed.returncode:
        detail = (completed.stderr or completed.stdout or "").strip()
        tail = "\n".join(detail.splitlines()[-20:])
        raise RuntimeError(
            f"modkit dmr pair failed (exit {completed.returncode}):\n{tail}"
        )
    return Path(kwargs.get("segment") or output)


def read_segments(path: str | Path) -> pd.DataFrame:
    """Parse a ``modkit dmr pair --segment`` file into a normalized frame.

    modkit 0.6 segmentation reports no p-value; significance comes from Cohen's
    h and its 95% CI (a segment is real when the CI excludes zero). Columns are
    matched loosely because names have shifted between modkit versions.
    """

    with open(path) as handle:
        first = handle.readline().rstrip("\n")
    if not first:
        return pd.DataFrame(
            columns=["chrom", "start", "end", "state", "num_sites", "effect",
                     "cohen_h", "cohen_h_low", "cohen_h_high",
                     "a_frac_modified", "b_frac_modified", "ci_excludes_zero"]
        )
    separator = "\t" if "\t" in first else r"\s+"
    frame = pd.read_csv(path, sep=separator, engine="python")
    frame.columns = [str(column).lstrip("#").strip().lower() for column in frame.columns]

    def pick(*names: str) -> str | None:
        for name in names:
            if name in frame.columns:
                return name
        return None

    c_chrom = pick("chrom", "chr", "chromosome")
    c_start = pick("chrom_start", "start", "region_start")
    c_end = pick("chrom_end", "end", "region_end")
    c_state = pick("name", "state")
    c_sites = pick("num_sites", "n_sites", "nsites")
    c_effect = pick("effect_size", "effect")
    result = pd.DataFrame({
        "chrom": frame[c_chrom].astype(str),
        "start": frame[c_start].astype(int),
        "end": frame[c_end].astype(int),
    })
    result["state"] = frame[c_state].astype(str).str.lower() if c_state else "different"
    result["num_sites"] = (
        pd.to_numeric(frame[c_sites], errors="coerce").fillna(0).astype(int) if c_sites else 0
    )
    result["effect"] = pd.to_numeric(frame[c_effect], errors="coerce") if c_effect else float("nan")
    for column in ("cohen_h", "cohen_h_low", "cohen_h_high", "a_frac_modified", "b_frac_modified"):
        result[column] = (
            pd.to_numeric(frame[column], errors="coerce") if column in frame.columns else float("nan")
        )
    if result["cohen_h_low"].notna().any():
        result["ci_excludes_zero"] = (
            (np.sign(result["cohen_h_low"]) == np.sign(result["cohen_h_high"]))
            & result["cohen_h_low"].notna()
        )
    else:
        result["ci_excludes_zero"] = True
    return result


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
