"""Read-level methylation extraction for candidate regions."""
from __future__ import annotations

import pysam
import pandas as pd


def methylation_tags(bam_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.has_tag("MM") or read.has_tag("Mm"):
                rows.append({
                    "read": read.query_name,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "haplotype": read.get_tag("HP") if read.has_tag("HP") else None,
                    "mm": read.get_tag("MM") if read.has_tag("MM") else read.get_tag("Mm"),
                })
    return pd.DataFrame(rows)
