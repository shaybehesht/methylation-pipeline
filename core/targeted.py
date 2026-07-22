"""Targeted candidate-gene methylation analysis (mirrors script 07).

For each promoter/gene-body region, take the CpGs covered at >= min_cov in all
three samples (so every comparison is paired on the same sites), then compute a
percentage-point effect size and a paired Wilcoxon signed-rank p-value. A region
is a proband-specific candidate when it differs from BOTH relatives by
>= min_delta pp with p < alpha, in the same direction, while the two relatives do
not differ there. chrX regions are compared to the female relative only.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import wilcoxon


def load_region(pileup_gz: str | Path, chrom: str, start: int, end: int, min_cov: int) -> dict[int, float]:
    """Return {position: percent_modified} for CpGs with valid coverage >= min_cov."""
    values: dict[int, float] = {}
    try:
        with pysam.TabixFile(str(pileup_gz)) as tabix:
            try:
                iterator = tabix.fetch(chrom, max(0, start), end)
            except ValueError:
                return values
            for line in iterator:
                fields = line.split("\t")
                if len(fields) < 11:
                    continue
                try:
                    if int(fields[9]) >= min_cov:
                        values[int(fields[1])] = float(fields[10])
                except ValueError:
                    continue
    except (OSError, ValueError):
        return values
    return values


def paired_compare(a_values, b_values) -> tuple[float, float]:
    """Return (delta in percentage points, Wilcoxon p-value) for paired vectors."""
    a = np.asarray(a_values, dtype=float)
    b = np.asarray(b_values, dtype=float)
    if len(a) < 3:
        return float("nan"), float("nan")
    delta = float(np.mean(a) - np.mean(b))
    difference = a - b
    if np.allclose(difference, 0):
        return delta, 1.0
    try:
        _, pvalue = wilcoxon(difference, zero_method="wilcox", alternative="two-sided")
    except ValueError:
        pvalue = float("nan")
    return delta, float(pvalue)


def score_regions(
    named_regions: pd.DataFrame,
    pileups: dict[str, Path],
    *,
    proband: str,
    relative_one: str,
    relative_two: str,
    female_relative: str | None = None,
    min_cov: int = 10,
    min_cpgs: int = 5,
    min_delta: float = 10.0,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Score every named region and flag proband-specific candidates.

    ``min_delta`` is in percentage points (matching script 07). Relatives are
    supplied by label; ``female_relative`` is the relative used for chrX
    comparisons (the proband's chrX is only comparable to a female relative).
    """

    rows: list[dict] = []
    for _, region in named_regions.iterrows():
        chrom, start, end = str(region["chrom"]), int(region["start"]), int(region["end"])
        per_sample = {
            label: load_region(pileups[label], chrom, start, end, min_cov)
            for label in (proband, relative_one, relative_two)
        }
        shared = sorted(
            set(per_sample[proband])
            & set(per_sample[relative_one])
            & set(per_sample[relative_two])
        )
        # Initialise the scoring columns to NaN so the frame always has them,
        # even when no region clears min_cpgs (otherwise the later column
        # references raise KeyError). NaN comparisons yield False -> no candidate.
        record = {
            "gene": region.get("gene", region.get("name", "")),
            "region": region.get("region", ""),
            "chrom": chrom, "start": start, "end": end,
            "n_cpgs": len(shared), "chrX": chrom == "chrX",
            "proband_meth": float("nan"),
            f"{relative_one}_meth": float("nan"),
            f"{relative_two}_meth": float("nan"),
            "delta_p_r1": float("nan"), "p_p_r1": float("nan"),
            "delta_p_r2": float("nan"), "p_p_r2": float("nan"),
            "delta_r1_r2": float("nan"), "p_r1_r2": float("nan"),
        }
        if len(shared) >= min_cpgs:
            proband_values = [per_sample[proband][pos] for pos in shared]
            r1_values = [per_sample[relative_one][pos] for pos in shared]
            r2_values = [per_sample[relative_two][pos] for pos in shared]
            record["proband_meth"] = float(np.mean(proband_values))
            record[f"{relative_one}_meth"] = float(np.mean(r1_values))
            record[f"{relative_two}_meth"] = float(np.mean(r2_values))
            record["delta_p_r1"], record["p_p_r1"] = paired_compare(proband_values, r1_values)
            record["delta_p_r2"], record["p_p_r2"] = paired_compare(proband_values, r2_values)
            record["delta_r1_r2"], record["p_r1_r2"] = paired_compare(r1_values, r2_values)
        rows.append(record)

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame

    enough = frame["n_cpgs"] >= min_cpgs
    sig_r1 = enough & (frame["delta_p_r1"].abs() >= min_delta) & (frame["p_p_r1"] < alpha)
    sig_r2 = enough & (frame["delta_p_r2"].abs() >= min_delta) & (frame["p_p_r2"] < alpha)
    concordant = np.sign(frame["delta_p_r1"]) == np.sign(frame["delta_p_r2"])
    null_quiet = ~((frame["delta_r1_r2"].abs() >= min_delta) & (frame["p_r1_r2"] < alpha))

    frame["candidate"] = sig_r1 & sig_r2 & concordant & null_quiet & ~frame["chrX"]
    # chrX: comparable to the female relative only.
    if female_relative == relative_one:
        frame.loc[frame["chrX"], "candidate"] = sig_r1[frame["chrX"]]
    elif female_relative == relative_two:
        frame.loc[frame["chrX"], "candidate"] = sig_r2[frame["chrX"]]

    frame["min_abs_delta"] = np.minimum(frame["delta_p_r1"].abs(), frame["delta_p_r2"].abs())
    if female_relative == relative_two:
        frame.loc[frame["chrX"], "min_abs_delta"] = frame.loc[frame["chrX"], "delta_p_r2"].abs()
    else:
        frame.loc[frame["chrX"], "min_abs_delta"] = frame.loc[frame["chrX"], "delta_p_r1"].abs()
    frame["direction"] = np.where(frame["delta_p_r1"] > 0, "hyper", "hypo")
    frame = frame.sort_values(["candidate", "min_abs_delta"], ascending=[False, False]).reset_index(drop=True)
    return frame
