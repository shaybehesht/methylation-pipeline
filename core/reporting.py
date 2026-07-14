"""Deterministic complete-run archive creation."""
from __future__ import annotations

import zipfile
from pathlib import Path

ARCHIVE_NAME = "methyl_trio_run.zip"


def build_run_archive(run_dir: str | Path, archive_name: str = ARCHIVE_NAME) -> Path:
    """Zip an entire run directory, excluding any previously built archive.

    Files are added in sorted order with fixed metadata so repeated builds over
    identical inputs produce byte-stable archives.
    """
    directory = Path(run_dir).expanduser().resolve()
    if not directory.is_dir():
        raise NotADirectoryError(directory)
    archive_path = directory / archive_name

    files = sorted(
        path for path in directory.rglob("*")
        if path.is_file() and path.name != archive_name
    )
    temporary = archive_path.with_suffix(".part")
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in files:
                info = zipfile.ZipInfo(str(path.relative_to(directory).as_posix()))
                info.date_time = (1980, 1, 1, 0, 0, 0)
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o644 << 16
                archive.writestr(info, path.read_bytes())
        temporary.replace(archive_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return archive_path
