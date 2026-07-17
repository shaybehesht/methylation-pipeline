from pathlib import Path

import pysam

from core.profile import parse_region, region_profile


def _write_pileup(path: Path, rows: list[tuple[int, float, int]]) -> Path:
    plain = path.with_suffix(".plain.bed")
    with plain.open("w") as handle:
        for pos, pct, cov in rows:
            handle.write("\t".join([
                "chr1", str(pos), str(pos + 1), "m", str(cov), "+",
                str(pos), str(pos + 1), "0,0,0", str(cov), f"{pct}",
                "0", "0", "0", "0", "0", "0", "0",
            ]) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)
    plain.unlink()
    return path


def test_parse_region_variants():
    assert parse_region("chr1:100-200") == ("chr1", 100, 200)
    assert parse_region("chr1:156,590,755-156,594,876") == ("chr1", 156590755, 156594876)
    assert parse_region("garbage") is None
    assert parse_region("chr1:200-100") is None


def test_region_profile_smooths_per_sample(tmp_path: Path):
    positions = list(range(1000, 1200, 10))
    a = _write_pileup(tmp_path / "P.bed.gz", [(p, 80.0, 20) for p in positions])
    b = _write_pileup(tmp_path / "M.bed.gz", [(p, 20.0, 20) for p in positions])
    profile = region_profile({"P": a, "M": b}, "chr1", 900, 1300, window=5)
    assert set(profile["sample"]) == {"P", "M"}
    assert {"pos", "pct", "smooth"} <= set(profile.columns)
    # smoothing keeps values within the observed range
    p_rows = profile[profile["sample"] == "P"]
    assert p_rows["smooth"].between(0, 100).all()
    assert len(p_rows) == len(positions)


def test_region_profile_empty_when_no_data(tmp_path: Path):
    empty = _write_pileup(tmp_path / "P.bed.gz", [])
    profile = region_profile({"P": empty}, "chr9", 1, 100, window=5)
    assert profile.empty
