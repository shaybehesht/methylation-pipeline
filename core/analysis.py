"""Trio-specific DMR intersection, ranking, and feasibility analysis."""
from __future__ import annotations

import numpy as np
import pandas as pd

REQUIRED = {"chrom", "start", "end", "effect", "pvalue", "n_sites"}


def _validate(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    missing = REQUIRED - set(frame.columns)
    if missing:
        raise ValueError(f"{name} is missing columns: {', '.join(sorted(missing))}")
    result = frame.copy()
    result["chrom"] = result["chrom"].astype(str)
    return result


def _overlaps(a: pd.Series, b: pd.DataFrame) -> pd.Series:
    return (
        (b["chrom"] == a["chrom"])
        & (b["start"] < a["end"])
        & (b["end"] > a["start"])
    )


def intersect_and_rank(
    proband_one: pd.DataFrame,
    proband_two: pd.DataFrame,
    relative_null: pd.DataFrame,
    *,
    null_percentile: float = 99,
    min_sites: int = 5,
    max_pval: float = 0.01,
    chromosome_validator=None,
) -> tuple[pd.DataFrame, float]:
    """Return concordant proband DMRs that exceed and do not overlap the null."""
    one = _validate(proband_one, "first proband comparison")
    two = _validate(proband_two, "second proband comparison")
    null = _validate(relative_null, "relative null comparison")
    null_effects = null["effect"].abs().dropna()
    cutoff = float(np.percentile(null_effects, null_percentile)) if len(null_effects) else 0.0
    null_variable = null[
        (null["effect"].abs() >= cutoff)
        & (null["n_sites"] >= min_sites)
        & (null["pvalue"] <= max_pval)
    ]

    rows: list[dict] = []
    for _, left in one.iterrows():
        if chromosome_validator and not chromosome_validator(str(left["chrom"])):
            continue
        matches = two.loc[_overlaps(left, two)]
        for _, right in matches.iterrows():
            if np.sign(left["effect"]) != np.sign(right["effect"]):
                continue
            if min(left["n_sites"], right["n_sites"]) < min_sites:
                continue
            if max(left["pvalue"], right["pvalue"]) > max_pval:
                continue
            if min(abs(left["effect"]), abs(right["effect"])) <= cutoff:
                continue
            start, end = max(int(left["start"]), int(right["start"])), min(int(left["end"]), int(right["end"]))
            candidate = pd.Series({"chrom": left["chrom"], "start": start, "end": end})
            if null_variable.loc[_overlaps(candidate, null_variable)].shape[0]:
                continue
            rows.append({
                "chrom": left["chrom"], "start": start, "end": end,
                "effect_1": float(left["effect"]), "effect_2": float(right["effect"]),
                "mean_abs_effect": float(np.mean([abs(left["effect"]), abs(right["effect"])])),
                "min_sites": int(min(left["n_sites"], right["n_sites"])),
                "max_pvalue": float(max(left["pvalue"], right["pvalue"])),
            })
    columns = ["chrom", "start", "end", "effect_1", "effect_2", "mean_abs_effect", "min_sites", "max_pvalue"]
    ranked = pd.DataFrame(rows, columns=columns)
    if not ranked.empty:
        ranked = ranked.sort_values(["mean_abs_effect", "max_pvalue"], ascending=[False, True]).reset_index(drop=True)
        ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))
    return ranked, cutoff


def feasibility(candidate_count: int, null_count: int) -> dict[str, float | str]:
    ratio = candidate_count / max(null_count, 1)
    if candidate_count >= 10 and ratio >= 2:
        verdict = "PURSUE"
    elif candidate_count >= 3 and ratio >= 0.75:
        verdict = "MARGINAL"
    else:
        verdict = "DO NOT PURSUE"
    return {"verdict": verdict, "ratio": ratio, "candidate_count": candidate_count, "null_count": null_count}


def targeted_test(values_a, values_b) -> dict[str, float]:
    from scipy.stats import mannwhitneyu

    a, b = np.asarray(values_a, dtype=float), np.asarray(values_b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if not len(a) or not len(b):
        return {"delta_pp": float("nan"), "pvalue": float("nan")}
    pvalue = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    return {"delta_pp": float((a.mean() - b.mean()) * 100), "pvalue": pvalue}
