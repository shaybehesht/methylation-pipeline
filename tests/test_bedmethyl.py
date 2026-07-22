from pathlib import Path
import os

import pytest

from core.bam_index import ensure_index
from core.bedmethyl import (
    bedmethyl_has_rows,
    pileup_modified_bases,
    validate_bedmethyl,
)
from core.dmr import run_pair
from core.pileup import build_command


def test_pileup_modified_bases_adds_5hmc_for_5mc():
    assert pileup_modified_bases(("5mC",)) == ("5mC", "5hmC")
    assert pileup_modified_bases(("5mC", "5hmC")) == ("5mC", "5hmC")
    assert pileup_modified_bases(("5hmC",)) == ("5hmC",)


def test_pileup_command_includes_both_mcs_when_configured_for_5mc_only():
    command = build_command("a.bam", "a.bed.gz", "ref.fa", modified_bases=pileup_modified_bases(("5mC",)))
    index = command.index("--modified-bases")
    assert command[index + 1:index + 3] == ["5mC", "5hmC"]


def test_validate_bedmethyl_detects_count_mismatch(tmp_path: Path):
    path = tmp_path / "pileup.bed"
    path.write_text(
        "\t".join([
            "chr14", "100", "101", "m", "10", "+", "100", "101", "0,0,0",
            "10", "60.0", "6", "3", "0", "0", "0", "0", "0",
        ]) + "\n",
        encoding="utf-8",
    )
    n_rows, n_invalid, examples = validate_bedmethyl(path)
    assert n_rows == 1
    assert n_invalid == 1
    assert "chr14:100" in examples[0]


def test_validate_bedmethyl_accepts_consistent_row(tmp_path: Path):
    path = tmp_path / "pileup.bed"
    path.write_text(
        "\t".join([
            "chr14", "100", "101", "m", "10", "+", "100", "101", "0,0,0",
            "10", "60.0", "6", "4", "0", "0", "0", "0", "0",
        ]) + "\n",
        encoding="utf-8",
    )
    n_rows, n_invalid, _ = validate_bedmethyl(path)
    assert n_rows == 1
    assert n_invalid == 0


def test_bedmethyl_has_rows(tmp_path: Path):
    empty = tmp_path / "empty.bed"
    empty.write_text("#comment\n", encoding="utf-8")
    assert bedmethyl_has_rows(empty) is False
    nonempty = tmp_path / "rows.bed"
    nonempty.write_text(
        "\t".join(["chr1", "1", "2", "m", "1", "+", "1", "2", "0", "1", "100", "1", "0", "0"]) + "\n",
        encoding="utf-8",
    )
    assert bedmethyl_has_rows(nonempty) is True


def test_run_pair_writes_empty_segments_for_sparse_chrom(tmp_path: Path):
    left = tmp_path / "left.bed.gz"
    right = tmp_path / "right.bed.gz"
    left.write_bytes(b"")
    right.write_bytes(b"")
    segment = tmp_path / "seg.bed"
    sites = tmp_path / "sites.bed"
    run_pair(str(left), str(right), str(sites), "ref.fa", segment=str(segment))
    assert segment.exists()
    assert "chrom" in segment.read_text(encoding="utf-8")


def test_run_pair_preflight_reports_invalid_rows(tmp_path: Path):
    left = tmp_path / "left.bed"
    right = tmp_path / "right.bed"
    row = "\t".join([
        "chr14", "100", "101", "m", "10", "+", "100", "101", "0,0,0",
        "10", "60.0", "6", "3", "0", "0", "0", "0", "0",
    ])
    left.write_text(row + "\n", encoding="utf-8")
    right.write_text(
        "\t".join([
            "chr14", "100", "101", "m", "10", "+", "100", "101", "0,0,0",
            "10", "60.0", "6", "4", "0", "0", "0", "0", "0",
        ]) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="invalid row"):
        run_pair(str(left), str(right), str(tmp_path / "out.bed"), "ref.fa")


def test_ensure_index_skips_existing_bai(tmp_path: Path, monkeypatch):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM")
    index = tmp_path / "sample.bam.bai"
    index.write_bytes(b"BAI")
    calls: list[list[str]] = []

    def fake_run(command, check):
        calls.append(command)

    monkeypatch.setattr("core.bam_index.subprocess.run", fake_run)
    ensure_index(str(bam))
    assert calls == []


def test_ensure_index_creates_missing_bai(tmp_path: Path, monkeypatch):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM")
    calls: list[list[str]] = []

    def fake_run(command, check):
        calls.append(command)
        Path(command[3]).write_bytes(b"BAI")

    monkeypatch.setattr("core.bam_index.subprocess.run", fake_run)
    ensure_index(str(bam))
    assert calls == [["samtools", "index", str(bam), str(tmp_path / "sample.bam.bai")]]


def test_ensure_index_skips_stale_symlink(tmp_path: Path, monkeypatch):
    bam = tmp_path / "sample.bam"
    target = tmp_path / "localized.bam"
    target.write_bytes(b"BAM2")
    bam.symlink_to(target)
    linked = tmp_path / "sample.bam.bai"
    old_bai = tmp_path / "localized.bam.bai"
    old_bai.write_bytes(b"OLD")
    linked.symlink_to(old_bai)
    os.utime(old_bai, (1, 1))
    os.utime(target, (2, 2))
    calls: list[list[str]] = []

    def fake_run(command, check):
        calls.append(command)

    monkeypatch.setattr("core.bam_index.subprocess.run", fake_run)
    ensure_index(str(bam))
    assert calls == []
