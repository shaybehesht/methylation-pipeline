from pathlib import Path

import pandas as pd

from core.annotations import extract_to_regions, panel_regions, write_bed3
from core.dmr import read_segments
from core.pileup import fold_hmc
from core.segments import (
    concordant_private,
    null_effect_threshold,
    proband_private_dmrs,
    significant,
    subtract_variable,
)


# ---------------------------------------------------------------- pileup fold
def test_fold_hmc_moves_n_other_into_canonical():
    # cols (1-indexed): 13=Ncanonical, 14=Nother
    row = "\t".join([
        "chr1", "100", "101", "m", "10", "+", "100", "101", "0,0,0",
        "10", "60.0", "6", "3", "1", "0", "0", "0", "0",
    ])
    folded = fold_hmc(row).split("\t")
    assert folded[12] == "4"   # 3 canonical + 1 other
    assert folded[13] == "0"   # other zeroed
    # valid coverage (col 10) and percent (col 11) untouched
    assert folded[9] == "10" and folded[10] == "60.0"


# ---------------------------------------------------------------- segments IO
SEGMENT_HEADER = (
    "chrom\tchrom_start\tchrom_end\tname\tscore\tnum_sites\teffect_size\t"
    "cohen_h\tcohen_h_low\tcohen_h_high\n"
)


def _write_segments(path: Path, rows: list[tuple]) -> Path:
    with path.open("w") as handle:
        handle.write(SEGMENT_HEADER)
        for row in rows:
            handle.write("\t".join(str(v) for v in row) + "\n")
    return path


def test_read_segments_parses_and_flags_ci(tmp_path: Path):
    path = _write_segments(tmp_path / "seg.bed", [
        ("chr1", 100, 200, "different", 10, 8, 0.40, 0.5, 0.3, 0.7),
        ("chr1", 300, 400, "same", 5, 6, 0.02, 0.1, -0.1, 0.2),
    ])
    frame = read_segments(path)
    assert list(frame["chrom"]) == ["chr1", "chr1"]
    assert frame.loc[0, "start"] == 100 and frame.loc[0, "end"] == 200
    assert frame.loc[0, "num_sites"] == 8
    assert bool(frame.loc[0, "ci_excludes_zero"]) is True    # 0.3..0.7 same sign
    assert bool(frame.loc[1, "ci_excludes_zero"]) is False   # -0.1..0.2 spans 0


def _frame(rows):
    return pd.DataFrame(
        rows, columns=["chrom", "start", "end", "state", "num_sites", "effect"]
    ).assign(ci_excludes_zero=True)


def test_null_threshold_uses_percentile_or_fallback():
    small = _frame([("chr1", 0, 1, "different", 8, 0.1)])
    assert null_effect_threshold(small, 99) == 0.25  # <50 segments -> fixed fallback
    big = _frame([("chr1", i, i + 1, "different", 8, 0.1 + i / 1000) for i in range(100)])
    assert null_effect_threshold(big, 90) > 0.1


def test_concordant_private_requires_overlap_and_direction():
    pm = _frame([("chr1", 100, 200, "different", 8, 0.40)])
    pb = _frame([("chr1", 120, 210, "different", 8, 0.35)])
    keep = concordant_private(pm, pb, effect=0.1, min_sites=5, require_ci=True)
    assert len(keep) == 1
    assert keep.loc[0, "start"] == 100 and keep.loc[0, "end"] == 210

    discordant = _frame([("chr1", 120, 210, "different", 8, -0.35)])
    assert concordant_private(pm, discordant, 0.1, 5, True).empty


def test_subtract_variable_removes_null_overlap():
    candidates = pd.DataFrame([{"chrom": "chr1", "start": 100, "end": 200}])
    null = _frame([("chr1", 150, 160, "different", 8, 0.5)])
    assert subtract_variable(candidates, null).empty


def test_proband_private_end_to_end_flags_imprinted():
    pm = _frame([
        ("chr1", 100, 200, "different", 8, 0.40),
        ("chr11", 1999500, 2001000, "different", 8, 0.45),  # H19 imprinted locus
    ])
    pb = _frame([
        ("chr1", 120, 210, "different", 8, 0.38),
        ("chr11", 1999600, 2001200, "different", 8, 0.44),
    ])
    mb = _frame([("chr2", 0, 10, "same", 8, 0.05)])
    ranked, cutoff = proband_private_dmrs(pm, pb, mb, null_percentile=90, min_sites=5)
    assert cutoff == 0.25  # <50 null segments -> fallback
    assert len(ranked) == 2
    assert set(ranked["chrom"]) == {"chr1", "chr11"}
    imprinted = ranked[ranked["chrom"] == "chr11"].iloc[0]["imprinted_control"]
    assert imprinted == "H19/ICR1"
    assert "rank" in ranked.columns and "mean_effect" in ranked.columns


def test_significant_requires_ci_when_asked():
    frame = _frame([("chr1", 0, 1, "different", 8, 0.5)])
    frame.loc[0, "ci_excludes_zero"] = False
    assert significant(frame, 0.1, 5, require_ci=True).empty
    assert not significant(frame, 0.1, 5, require_ci=False).empty


# ---------------------------------------------------------------- annotations
GTF = (
    '#comment\n'
    'chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tgene_name "HESX1";\n'
    'chr3\tsrc\tgene\t2000\t8000\t.\t-\t.\tgene_name "SOX2";\n'
)


def test_panel_regions_builds_promoter_body_and_extract(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF, encoding="utf-8")
    named, extract, missing = panel_regions(gtf, ["HESX1", "SOX2", "NOPE"], 2000, 5000)
    assert missing == ["NOPE"]
    assert set(named["region"]) == {"promoter", "body"}
    assert set(named["gene"]) == {"HESX1", "SOX2"}
    # promoter for + strand HESX1 centres on TSS = start (999, 0-based)
    hesx1_prom = named[(named["gene"] == "HESX1") & (named["region"] == "promoter")].iloc[0]
    assert hesx1_prom["start"] == max(0, 999 - 2000)
    assert list(extract.columns) == ["chrom", "start", "end"]


def test_write_bed3_emits_three_columns(tmp_path: Path):
    frame = pd.DataFrame([{"chrom": "chr1", "start": 10, "end": 20, "name": "x"}])
    path = write_bed3(frame, tmp_path / "scope.bed")
    lines = path.read_text().strip().splitlines()
    assert lines == ["chr1\t10\t20"]


def test_extract_to_regions_converts_bed_to_1based():
    frame = pd.DataFrame([
        {"chrom": "chr3", "start": 100, "end": 200},
        {"chrom": "chrX", "start": 0, "end": 50},
    ])
    assert extract_to_regions(frame) == ["chr3:101-200", "chrX:1-50"]
