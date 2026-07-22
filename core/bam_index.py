"""Ensure modBAM indexes are present and newer than their BAM."""
from __future__ import annotations

import subprocess
from pathlib import Path


def ensure_index(bam_path: str, log=None) -> None:
    """Create or refresh ``.bai`` when missing or older than the BAM."""
    bam = Path(bam_path)
    if not bam.exists():
        raise FileNotFoundError(bam_path)
    index = Path(f"{bam}.bai")
    needs_index = not index.exists()
    if not needs_index and index.stat().st_mtime < bam.stat().st_mtime:
        needs_index = True
        if log is not None:
            log.write(f"Re-indexing stale BAI (older than BAM): {bam}\n")
    if needs_index:
        if log is not None and index.exists():
            pass  # message already written for stale case
        elif log is not None:
            log.write(f"Indexing BAM: {bam}\n")
        subprocess.run(["samtools", "index", str(bam)], check=True)
