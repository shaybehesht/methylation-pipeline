"""Ensure modBAM indexes are present and newer than their BAM."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _writable_index_path(bam: Path) -> Path:
    """Return a local ``.bai`` path Cromwell/Terra can write."""
    index = Path(f"{bam}.bai")
    if index.exists() and not index.is_symlink() and os.access(index, os.W_OK):
        return index
    # Localized BAMs are often symlinks with a sibling symlinked BAI on read-only
    # storage; write (or refresh) an index next to the staged BAM instead.
    return bam.with_name(bam.name + ".bai")


def ensure_index(bam_path: str, log=None) -> None:
    """Create a missing ``.bai``, or refresh one we can write locally.

    Stale indexes (older than the BAM) only emit a warning. Re-indexing multi-GB
    PacBio modBAMs on Terra is slow, memory-heavy, and often impossible when the
    localized BAM lives on a read-only mount — htslib still reads the BAM fine.
    """
    bam = Path(bam_path)
    if not bam.exists():
        raise FileNotFoundError(bam_path)
    index = _writable_index_path(bam)
    linked = Path(f"{bam}.bai")
    if linked.exists() and linked.is_symlink() and linked.resolve() != index.resolve():
        if log is not None:
            log.write(f"Using writable BAI path {index} (localized index is a symlink)\n")
    if linked.exists() and linked.stat().st_mtime < bam.stat().st_mtime:
        if log is not None:
            log.write(
                f"Warning: BAI older than BAM for {bam} "
                f"(continuing without re-index; refresh the index in your bucket if needed)\n"
            )
    if index.exists():
        return
    if log is not None:
        log.write(f"Indexing BAM -> {index}\n")
    if linked.exists() and linked.is_symlink():
        linked.unlink()
    index.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["samtools", "index", str(bam), str(index)], check=True)
