"""Ensure modBAM indexes are present and newer than their BAM."""
from __future__ import annotations

import subprocess
from pathlib import Path


def ensure_index(bam_path: str, log=None) -> None:
    """Use an existing ``.bai`` when provided; only build one if it is missing.

    WDL/Terra runs always pass ``*_bai`` alongside each BAM. Re-indexing large
    PacBio modBAMs on read-only localized storage is slow and unnecessary — the
    htslib stale-index warning is harmless for pileup.
    """
    bam = Path(bam_path)
    if not bam.exists():
        raise FileNotFoundError(bam_path)
    index = Path(f"{bam}.bai")
    if index.exists():
        if index.stat().st_mtime < bam.stat().st_mtime and log is not None:
            log.write(
                f"Note: BAI older than BAM for {bam.name} "
                f"(using your supplied index; re-index in the bucket only if reads look wrong)\n"
            )
        return
    if log is not None:
        log.write(f"No BAI found beside {bam}; running samtools index\n")
    subprocess.run(["samtools", "index", str(bam), str(index)], check=True)
