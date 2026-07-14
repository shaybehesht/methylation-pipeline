import os
from pathlib import Path

import pytest

from app import file_picker


@pytest.fixture
def rooted(tmp_path, monkeypatch):
    monkeypatch.setenv("METHYL_TRIO_DATA_ROOT", str(tmp_path))
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "sample.bam").write_bytes(b"data")
    return tmp_path


def test_data_root_uses_environment(rooted):
    assert file_picker.data_root() == rooted.resolve()


def test_resolve_in_root_accepts_child(rooted):
    resolved = file_picker.resolve_in_root("sub/sample.bam")
    assert resolved == (rooted / "sub" / "sample.bam").resolve()


def test_resolve_in_root_blocks_traversal(rooted):
    assert file_picker.resolve_in_root("../etc/passwd") is None
    assert file_picker.resolve_in_root("missing.bam") is None


def test_bam_index_detects_both_conventions(tmp_path):
    bam = tmp_path / "a.bam"
    bam.write_bytes(b"x")
    assert file_picker.bam_index(str(bam)) is None
    companion = tmp_path / "a.bam.bai"
    companion.write_bytes(b"i")
    assert file_picker.bam_index(str(bam)) == str(companion)
    companion.unlink()
    (tmp_path / "a.bai").write_bytes(b"i")
    assert file_picker.bam_index(str(bam)).endswith("a.bai")
