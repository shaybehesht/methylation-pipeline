"""Offline GTF gene lookup and target BED generation."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

GENE_NAME = re.compile(r'gene_name "([^"]+)"')


def genes_from_gtf(path: str | Path, symbols: list[str], promoter_pad: int, body_pad: int) -> pd.DataFrame:
    wanted = {symbol.strip().upper() for symbol in symbols if symbol.strip()}
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = GENE_NAME.search(fields[8])
            if not match or match.group(1).upper() not in wanted:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            tss = start if strand == "+" else end
            rows.append({
                "chrom": fields[0],
                "start": max(0, min(start - body_pad, tss - promoter_pad)),
                "end": max(end + body_pad, tss + promoter_pad),
                "gene": match.group(1),
                "strand": strand,
            })
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "strand"])


def write_bed(frame: pd.DataFrame, path: str | Path) -> Path:
    output = Path(path)
    frame[["chrom", "start", "end", "gene"]].to_csv(output, sep="\t", header=False, index=False)
    return output
