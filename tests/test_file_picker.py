import os
from pathlib import Path

from app import file_picker
from app.file_picker import (
    BAM_EXTENSIONS,
    data_root,
    data_roots,
    detect_bam_index,
    list_entries,
    resolve_within,
)


def test_resolve_within_confines_traversal(tmp_path: Path):
    root = tmp_path / "data"
    (root / "sub").mkdir(parents=True)
    inside = root / "sub"

    assert resolve_within(root, inside) == inside.resolve()
    assert resolve_within(root, root) == root.resolve()
    # Classic escape attempts must be rejected.
    assert resolve_within(root, root / ".." / "secret") is None
    assert resolve_within(root, tmp_path / "outside") is None


def test_resolve_within_rejects_symlink_escape(tmp_path: Path):
    root = tmp_path / "data"
    root.mkdir()
    secret = tmp_path / "secret.txt"
    secret.write_text("classified", encoding="utf-8")
    link = root / "escape"
    link.symlink_to(secret)

    # The symlink target resolves outside the root and must be refused.
    assert resolve_within(root, link) is None


def test_list_entries_filters_extensions_and_hidden(tmp_path: Path):
    (tmp_path / "a.bam").write_bytes(b"")
    (tmp_path / "b.txt").write_bytes(b"")
    (tmp_path / ".hidden.bam").write_bytes(b"")
    (tmp_path / "nested").mkdir()
    (tmp_path / ".secret").mkdir()

    dirs, files = list_entries(tmp_path, BAM_EXTENSIONS)
    assert [d.name for d in dirs] == ["nested"]
    assert [f.name for f in files] == ["a.bam"]


def test_detect_bam_index_variants(tmp_path: Path):
    dotbai = tmp_path / "sample_a.bam"
    dotbai.write_bytes(b"")
    (tmp_path / "sample_a.bam.bai").write_bytes(b"")
    assert detect_bam_index(dotbai) == tmp_path / "sample_a.bam.bai"

    stripped = tmp_path / "sample_b.bam"
    stripped.write_bytes(b"")
    (tmp_path / "sample_b.bai").write_bytes(b"")
    assert detect_bam_index(stripped) == tmp_path / "sample_b.bai"

    missing = tmp_path / "sample_c.bam"
    missing.write_bytes(b"")
    assert detect_bam_index(missing) is None


def test_data_root_uses_env(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(file_picker, "_EXTERNAL_MOUNT_PARENTS", ())
    monkeypatch.setenv("METHYL_TRIO_DATA_ROOT", str(tmp_path))
    assert data_root() == tmp_path.resolve()

    monkeypatch.delenv("METHYL_TRIO_DATA_ROOT", raising=False)
    assert data_root() == Path.home().resolve()


def test_data_roots_accepts_multiple_locations(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(file_picker, "_EXTERNAL_MOUNT_PARENTS", ())
    first = tmp_path / "local"
    second = tmp_path / "drive"
    first.mkdir()
    second.mkdir()
    monkeypatch.setenv("METHYL_TRIO_DATA_ROOT", os.pathsep.join([str(first), str(second)]))
    assert data_roots() == [first.resolve(), second.resolve()]


def test_data_roots_auto_detect_external_drives(monkeypatch, tmp_path: Path):
    mount_parent = tmp_path / "Volumes"
    (mount_parent / "MyDrive").mkdir(parents=True)
    empty_parent = tmp_path / "media"
    empty_parent.mkdir()
    monkeypatch.setattr(
        file_picker, "_EXTERNAL_MOUNT_PARENTS", (str(mount_parent), str(empty_parent))
    )
    monkeypatch.delenv("METHYL_TRIO_DATA_ROOT", raising=False)

    roots = data_roots()
    assert roots[0] == Path.home().resolve()
    # The populated mount parent is offered; the empty one is skipped.
    assert mount_parent.resolve() in roots
    assert empty_parent.resolve() not in roots
