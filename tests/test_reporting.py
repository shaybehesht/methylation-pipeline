import zipfile
from pathlib import Path

import pytest

from core.reporting import archive_size, build_run_archive, human_size


def _make_run(tmp_path: Path) -> Path:
    run_dir = tmp_path / "run-1"
    (run_dir / "figures").mkdir(parents=True)
    (run_dir / "tables").mkdir()
    (run_dir / "config.json").write_text('{"assembly": "hg38"}', encoding="utf-8")
    (run_dir / "pipeline.log").write_text("done\n", encoding="utf-8")
    (run_dir / "Proband.bed.gz").write_bytes(b"pileup")
    (run_dir / "tables" / "proband_specific_DMRs.tsv").write_text("rank\n1\n", encoding="utf-8")
    (run_dir / "figures" / "dmr_effects.png").write_bytes(b"\x89PNG\r\n")
    (run_dir / "report.html").write_text("<html></html>", encoding="utf-8")
    return run_dir


def test_archive_contains_full_run_and_excludes_itself(tmp_path: Path):
    run_dir = _make_run(tmp_path)
    archive = build_run_archive(run_dir)

    assert archive == run_dir / "complete_run.zip"
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())

    assert names == {
        "config.json",
        "pipeline.log",
        "Proband.bed.gz",
        "tables/proband_specific_DMRs.tsv",
        "figures/dmr_effects.png",
        "report.html",
    }
    assert "complete_run.zip" not in names


def test_archive_is_deterministic(tmp_path: Path):
    run_dir = _make_run(tmp_path)
    first = build_run_archive(run_dir, tmp_path / "a.zip")
    second = build_run_archive(run_dir, tmp_path / "b.zip")
    assert first.read_bytes() == second.read_bytes()


def test_rebuild_excludes_prior_archive(tmp_path: Path):
    run_dir = _make_run(tmp_path)
    build_run_archive(run_dir)
    # Rebuilding must not embed the previous archive into the new one.
    archive = build_run_archive(run_dir)
    with zipfile.ZipFile(archive) as bundle:
        assert "complete_run.zip" not in bundle.namelist()


def test_missing_run_directory_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        build_run_archive(tmp_path / "missing")


def test_size_helpers(tmp_path: Path):
    run_dir = _make_run(tmp_path)
    archive = build_run_archive(run_dir)
    assert archive_size(archive) > 0
    assert archive_size(tmp_path / "nope.zip") == 0
    assert human_size(0) == "0 B"
    assert human_size(1536).endswith("KB")
    assert human_size(5 * 1024 * 1024).endswith("MB")
