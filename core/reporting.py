"""Deterministic packaging of a complete run directory into a single ZIP."""
from __future__ import annotations

import os
import zipfile
from pathlib import Path

DEFAULT_ARCHIVE_NAME = "complete_run.zip"
# ZIP timestamps only reach back to 1980; use it so identical inputs produce
# byte-identical archives regardless of when they are built.
_FIXED_TIMESTAMP = (1980, 1, 1, 0, 0, 0)


def _iter_run_files(run_dir: Path, archive_path: Path) -> list[Path]:
    archive_resolved = archive_path.resolve()
    files = []
    for path in run_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.resolve() == archive_resolved:
            continue
        if path.name.endswith(".part"):
            continue
        files.append(path)
    return sorted(files, key=lambda item: item.relative_to(run_dir).as_posix())


def build_run_archive(
    run_dir: str | Path, archive_path: str | Path | None = None
) -> Path:
    """ZIP the entire ``run_dir`` deterministically, excluding the archive itself.

    Every regular file under ``run_dir`` (manifests, logs, pileups, pairwise
    outputs, tables, figures, and the HTML report) is included with a stable
    ordering, timestamp, and permission bits so repeated builds of an unchanged
    run yield identical bytes. The archive is written atomically.
    """

    run_path = Path(run_dir)
    if not run_path.is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {run_path}")
    archive = Path(archive_path) if archive_path else run_path / DEFAULT_ARCHIVE_NAME

    files = _iter_run_files(run_path, archive)
    partial = archive.with_name(archive.name + ".part")
    try:
        with zipfile.ZipFile(partial, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
            for file_path in files:
                arcname = file_path.relative_to(run_path).as_posix()
                info = zipfile.ZipInfo(arcname, date_time=_FIXED_TIMESTAMP)
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o644 << 16
                bundle.writestr(info, file_path.read_bytes())
        os.replace(partial, archive)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    return archive


def archive_size(path: str | Path) -> int:
    """Return the archive size in bytes (0 when it does not exist)."""

    archive = Path(path)
    return archive.stat().st_size if archive.exists() else 0


def human_size(num_bytes: int) -> str:
    """Format a byte count as a short human-readable string."""

    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"
