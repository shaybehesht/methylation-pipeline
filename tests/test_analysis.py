import pandas as pd

from core.analysis import feasibility, intersect_and_rank, targeted_test


def frame(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "effect", "pvalue", "n_sites"])


def test_intersection_requires_concordance_cutoff_and_null_subtraction():
    one = frame([
        ("chr1", 100, 200, 0.40, 0.001, 8),
        ("chr1", 300, 400, -0.35, 0.001, 8),
        ("chrX", 500, 600, 0.50, 0.001, 8),
    ])
    two = frame([
        ("chr1", 120, 210, 0.38, 0.002, 7),
        ("chr1", 320, 410, 0.37, 0.001, 8),
        ("chrX", 510, 610, 0.48, 0.001, 8),
    ])
    null = frame([
        ("chr2", 0, 10, 0.10, 0.5, 8),
        ("chr1", 130, 140, 0.05, 0.5, 8),
    ])
    result, cutoff = intersect_and_rank(
        one, two, null, null_percentile=90,
        chromosome_validator=lambda chrom: chrom != "chrX",
    )
    assert cutoff == 0.095
    assert result.empty  # first overlaps null; second is discordant; chrX invalid


def test_ranked_candidate_and_feasibility():
    one = frame([("chr1", 100, 200, 0.40, 0.001, 8)])
    two = frame([("chr1", 120, 210, 0.30, 0.002, 7)])
    null = frame([("chr2", 0, 10, 0.05, 0.5, 8)])
    result, _ = intersect_and_rank(one, two, null)
    assert result.loc[0, "rank"] == 1
    assert result.loc[0, "start"] == 120
    assert feasibility(12, 4)["verdict"] == "PURSUE"
    assert feasibility(1, 10)["verdict"] == "DO NOT PURSUE"


def test_targeted_wilcoxon_summary():
    result = targeted_test([0.8, 0.9, 0.85], [0.1, 0.2, 0.15])
    assert result["delta_pp"] == 70
    assert 0 <= result["pvalue"] <= 1
