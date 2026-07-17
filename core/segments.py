"""Genome-wide segmentation analysis (mirrors scripts 02 and 11).

Given the three ``modkit dmr pair --segment`` outputs (proband-vs-relative-1,
proband-vs-relative-2, relative-1-vs-relative-2), find regions where the proband
differs from BOTH relatives in the same direction, using the relative-vs-relative
comparison as the empirical null, and subtracting regions that are merely
variable between the two relatives. Significance uses effect size, site count,
and the Cohen's h 95% CI (which must exclude zero), never a p-value.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Known imprinted gDMRs (GRCh38). These should sit near 50% methylation and be
# roughly equal across a trio; if they top the ranked list, suspect the pipeline
# rather than the biology.
IMPRINTED_GRCH38 = [
    ("chr11", 1999000, 2003000, "H19/ICR1"),
    ("chr11", 2698000, 2701000, "KCNQ1OT1/ICR2"),
    ("chr15", 24954000, 24956000, "SNRPN"),
    ("chr15", 23565000, 23575000, "MKRN3/MAGEL2"),
    ("chr7", 130490000, 130495000, "MEST"),
    ("chr20", 58838000, 58851000, "GNAS"),
    ("chr14", 100810000, 100830000, "MEG3/DLK1"),
    ("chr6", 143990000, 144010000, "PLAGL1"),
    ("chr19", 53190000, 53200000, "PEG3"),
    ("chr16", 3480000, 3495000, "ZNF597"),
]


def flip_segments(segments: pd.DataFrame) -> pd.DataFrame:
    """Turn an 'A vs B' segment frame into 'B vs A' (negate effect and Cohen's h CI)."""
    if segments.empty:
        return segments.copy()
    flipped = segments.copy()
    flipped["effect"] = -flipped["effect"]
    for column in ("cohen_h",):
        if column in flipped.columns:
            flipped[column] = -flipped[column]
    if "cohen_h_low" in flipped.columns and "cohen_h_high" in flipped.columns:
        low, high = flipped["cohen_h_low"].copy(), flipped["cohen_h_high"].copy()
        flipped["cohen_h_low"], flipped["cohen_h_high"] = -high, -low
    if "a_frac_modified" in flipped.columns and "b_frac_modified" in flipped.columns:
        flipped[["a_frac_modified", "b_frac_modified"]] = flipped[
            ["b_frac_modified", "a_frac_modified"]
        ].to_numpy()
    return flipped


def private_count(
    proband_one: pd.DataFrame, proband_two: pd.DataFrame, relative_null: pd.DataFrame,
    effect: float, min_sites: int, require_ci: bool = True,
) -> int:
    """Count private DMRs at a fixed effect threshold (used by the sweep)."""
    candidates = concordant_private(proband_one, proband_two, effect, min_sites, require_ci)
    candidates = subtract_variable(
        candidates, significant(relative_null, effect, min_sites, require_ci)
    )
    return len(candidates)


def significant(segments: pd.DataFrame, effect: float, min_sites: int, require_ci: bool) -> pd.DataFrame:
    if segments.empty:
        return segments
    mask = (
        segments["state"].str.contains("diff", na=False)
        & (segments["num_sites"] >= min_sites)
        & (segments["effect"].abs() >= effect)
    )
    if require_ci and "ci_excludes_zero" in segments.columns:
        mask &= segments["ci_excludes_zero"]
    return segments[mask]


def null_effect_threshold(null_segments: pd.DataFrame, percentile: float) -> float:
    """Return the given percentile of |effect| across the relative-null segments."""
    effects = null_segments["effect"].abs().dropna() if not null_segments.empty else pd.Series(dtype=float)
    if len(effects) < 50:
        return 0.25
    return float(np.percentile(effects, percentile))


def _overlaps(row: pd.Series, frame: pd.DataFrame) -> pd.Series:
    return (
        (frame["chrom"] == row["chrom"])
        & (frame["start"] < row["end"])
        & (frame["end"] > row["start"])
    )


def concordant_private(
    proband_one: pd.DataFrame,
    proband_two: pd.DataFrame,
    effect: float,
    min_sites: int,
    require_ci: bool,
    chromosome_validator=None,
) -> pd.DataFrame:
    """Regions significant in both proband comparisons with a concordant sign."""
    left = significant(proband_one, effect, min_sites, require_ci)
    right = significant(proband_two, effect, min_sites, require_ci)
    if left.empty or right.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "effect_1", "effect_2",
                                     "mean_effect", "min_abs_effect", "num_sites"])
    right_by_chrom = {chrom: group for chrom, group in right.groupby("chrom")}
    rows: list[dict] = []
    for _, region in left.iterrows():
        if chromosome_validator and not chromosome_validator(str(region["chrom"])):
            continue
        candidates = right_by_chrom.get(region["chrom"])
        if candidates is None:
            continue
        hit = candidates[
            (candidates["start"] < region["end"])
            & (candidates["end"] > region["start"])
            & (np.sign(candidates["effect"]) == np.sign(region["effect"]))
        ]
        if hit.empty:
            continue
        rows.append({
            "chrom": region["chrom"],
            "start": int(min(region["start"], hit["start"].min())),
            "end": int(max(region["end"], hit["end"].max())),
            "effect_1": float(region["effect"]),
            "effect_2": float(hit["effect"].mean()),
            "mean_effect": float((region["effect"] + hit["effect"].mean()) / 2),
            "min_abs_effect": float(min(abs(region["effect"]), hit["effect"].abs().min())),
            "num_sites": int(max(region["num_sites"], hit["num_sites"].max())),
        })
    if not rows:
        return pd.DataFrame(columns=["chrom", "start", "end", "effect_1", "effect_2",
                                     "mean_effect", "min_abs_effect", "num_sites"])
    return pd.DataFrame(rows).drop_duplicates(subset=["chrom", "start", "end"])


def subtract_variable(candidates: pd.DataFrame, null_variable: pd.DataFrame) -> pd.DataFrame:
    """Drop candidates that overlap a region also variable between the relatives."""
    if candidates.empty or null_variable.empty:
        return candidates
    keep = []
    for index, region in candidates.iterrows():
        if not _overlaps(region, null_variable).any():
            keep.append(index)
    return candidates.loc[keep]


def annotate_imprinted(candidates: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        result = candidates.copy()
        result["imprinted_control"] = pd.Series(dtype=str)
        return result
    labels = []
    for _, region in candidates.iterrows():
        label = ""
        for chrom, start, end, name in IMPRINTED_GRCH38:
            if region["chrom"] == chrom and region["start"] < end and region["end"] > start:
                label = name
                break
        labels.append(label)
    result = candidates.copy()
    result["imprinted_control"] = labels
    return result


def proband_private_dmrs(
    proband_one: pd.DataFrame,
    proband_two: pd.DataFrame,
    relative_null: pd.DataFrame,
    *,
    null_percentile: float = 99,
    min_sites: int = 5,
    require_ci: bool = True,
    chromosome_validator=None,
) -> tuple[pd.DataFrame, float]:
    """Return ranked proband-private DMRs and the empirical-null effect cutoff."""
    cutoff = null_effect_threshold(relative_null, null_percentile)
    candidates = concordant_private(
        proband_one, proband_two, cutoff, min_sites, require_ci, chromosome_validator
    )
    null_variable = significant(relative_null, cutoff, min_sites, require_ci)
    candidates = subtract_variable(candidates, null_variable)
    candidates = annotate_imprinted(candidates)
    if not candidates.empty:
        candidates["direction"] = np.where(candidates["mean_effect"] > 0, "hyper", "hypo")
        candidates["length"] = candidates["end"] - candidates["start"]
        candidates = candidates.sort_values("min_abs_effect", ascending=False).reset_index(drop=True)
        candidates.insert(0, "rank", np.arange(1, len(candidates) + 1))
    return candidates, cutoff
