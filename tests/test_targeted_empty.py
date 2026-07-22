"""Regression tests for targeted scoring when regions lack sufficient coverage."""
from pathlib import Path

import pandas as pd

import core.targeted as targeted


def _regions():
    return pd.DataFrame([
        {"gene": "SNRPN", "region": "promoter", "chrom": "chr15", "start": 100, "end": 300},
        {"gene": "MECP2", "region": "promoter", "chrom": "chrX", "start": 100, "end": 300},
    ])


def test_score_regions_no_shared_cpgs_does_not_crash(monkeypatch):
    # Simulate no covered CpGs anywhere (the PacBio failure mode): every region
    # falls below min_cpgs, so the scoring columns were previously never created.
    monkeypatch.setattr(targeted, "load_region", lambda *a, **k: {})
    scored = targeted.score_regions(
        _regions(),
        {"P": Path("p.bed.gz"), "R1": Path("r1.bed.gz"), "R2": Path("r2.bed.gz")},
        proband="P", relative_one="R1", relative_two="R2",
    )
    assert not scored.empty
    # columns must exist even though nothing cleared the threshold
    for col in ("delta_p_r1", "p_p_r1", "delta_p_r2", "delta_r1_r2", "candidate"):
        assert col in scored.columns
    assert int(scored["candidate"].sum()) == 0


def test_score_regions_with_coverage_flags_candidate(monkeypatch):
    def fake_load(pileup, chrom, start, end, min_cov):
        name = str(pileup)
        # proband hypermethylated vs both relatives at 8 shared CpGs
        base = {"p.bed.gz": 90.0, "r1.bed.gz": 8.0, "r2.bed.gz": 7.0}[name]
        return {pos: base + (pos % 3) for pos in range(start, start + 8)}

    monkeypatch.setattr(targeted, "load_region", fake_load)
    scored = targeted.score_regions(
        _regions(),
        {"P": Path("p.bed.gz"), "R1": Path("r1.bed.gz"), "R2": Path("r2.bed.gz")},
        proband="P", relative_one="R1", relative_two="R2",
        female_relative="R1", min_cov=5, min_cpgs=5, min_delta=10.0, alpha=0.05,
    )
    assert "delta_p_r1" in scored.columns
    # the autosomal SNRPN region should be a strong proband-hyper candidate
    snrpn = scored[scored["gene"] == "SNRPN"].iloc[0]
    assert snrpn["delta_p_r1"] > 10
    assert bool(scored["candidate"].any())
