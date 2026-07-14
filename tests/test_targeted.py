from pathlib import Path

import pandas as pd
import pysam

from core.targeted import paired_compare, score_regions


def _write_pileup(path: Path, positions: list[tuple[int, float, int]]) -> Path:
    """Write a minimal modkit-style bedMethyl and bgzip+tabix it.

    Each row: chrom start end mod score strand start end color coverage(10) percent(11)...
    """
    plain = path.with_suffix(".plain.bed")
    with plain.open("w") as handle:
        for pos, pct, cov in positions:
            handle.write("\t".join([
                "chr1", str(pos), str(pos + 1), "m", str(cov), "+",
                str(pos), str(pos + 1), "0,0,0", str(cov), f"{pct}",
                "0", "0", "0", "0", "0", "0", "0",
            ]) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)
    plain.unlink()
    return path


def test_paired_compare_delta_and_pvalue():
    delta, pvalue = paired_compare([80.0, 90.0, 85.0, 88.0], [10.0, 20.0, 15.0, 12.0])
    assert round(delta, 1) == 71.5
    assert 0 <= pvalue <= 1


def test_score_regions_flags_proband_specific_candidate(tmp_path: Path):
    positions = [100, 130, 160, 190, 220, 250]
    proband = _write_pileup(tmp_path / "P.bed.gz", [(p, 85.0, 20) for p in positions])
    r1 = _write_pileup(tmp_path / "M.bed.gz", [(p, 20.0, 20) for p in positions])
    r2 = _write_pileup(tmp_path / "B.bed.gz", [(p, 22.0, 20) for p in positions])

    regions = pd.DataFrame([
        {"chrom": "chr1", "start": 90, "end": 300, "name": "GENEA|promoter",
         "gene": "GENEA", "region": "promoter"},
    ])
    scored = score_regions(
        regions, {"P": proband, "M": r1, "B": r2},
        proband="P", relative_one="M", relative_two="B", female_relative="M",
        min_cov=10, min_cpgs=5, min_delta=10.0, alpha=0.05,
    )
    row = scored.iloc[0]
    assert row["n_cpgs"] == len(positions)
    assert bool(row["candidate"]) is True
    assert row["direction"] == "hyper"
    assert abs(row["delta_p_r1"] - 65.0) < 1e-6


def test_score_regions_not_candidate_when_proband_tracks_relatives(tmp_path: Path):
    positions = [100, 130, 160, 190, 220, 250]
    proband = _write_pileup(tmp_path / "P.bed.gz", [(p, 21.0, 20) for p in positions])
    r1 = _write_pileup(tmp_path / "M.bed.gz", [(p, 20.0, 20) for p in positions])
    r2 = _write_pileup(tmp_path / "B.bed.gz", [(p, 22.0, 20) for p in positions])
    regions = pd.DataFrame([
        {"chrom": "chr1", "start": 90, "end": 300, "name": "GENEA|body",
         "gene": "GENEA", "region": "body"},
    ])
    scored = score_regions(
        regions, {"P": proband, "M": r1, "B": r2},
        proband="P", relative_one="M", relative_two="B", female_relative="M",
    )
    assert bool(scored.iloc[0]["candidate"]) is False
