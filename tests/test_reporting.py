import zipfile
from pathlib import Path

import pytest

from core.reporting import build_run_archive


def _populate(run_dir: Path) -> None:
    (run_dir / "config.json").write_text("{}")
    (run_dir / "pipeline.log").write_text("log")
    (run_dir / "proband_specific_DMRs.tsv").write_text("rank\n1\n")
    (run_dir / "report.html").write_text("<html></html>")
    (run_dir / "sub").mkdir()
    (run_dir / "sub" / "Proband.bed.gz").write_bytes(b"pileup")


def test_archive_includes_everything_and_excludes_itself(tmp_path):
    _populate(tmp_path)
    archive = build_run_archive(tmp_path)
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
    assert {
        "config.json", "pipeline.log", "proband_specific_DMRs.tsv",
        "report.html", "sub/Proband.bed.gz",
    } <= names
    assert "methyl_trio_run.zip" not in names


def test_archive_is_deterministic(tmp_path):
    _populate(tmp_path)
    first = build_run_archive(tmp_path).read_bytes()
    second = build_run_archive(tmp_path).read_bytes()
    assert first == second


def test_archive_requires_directory(tmp_path):
    missing = tmp_path / "nope"
    with pytest.raises(NotADirectoryError):
        build_run_archive(missing)
