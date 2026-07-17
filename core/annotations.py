"""Offline GTF gene lookup and target BED generation."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

GENE_NAME = re.compile(r'gene_name "([^"]+)"')

# GENCODE renamed some symbols; alias current names back to requested ones.
GENE_ALIASES = {"RETREG1": "FAM134B"}


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


def all_genes(gtf: str | Path) -> pd.DataFrame:
    """Load every gene (chrom, start, end, tss, gene) from a GENCODE GTF."""
    rows = []
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = GENE_NAME.search(fields[8])
            if not match:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            tss = start if strand == "+" else end
            rows.append({"chrom": fields[0], "start": start, "end": end, "tss": tss, "gene": match.group(1)})
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "tss", "gene"])


def annotate_with_genes(
    candidates: pd.DataFrame, genes: pd.DataFrame, promoter_pad: int = 2000, limit: int = 6
) -> pd.DataFrame:
    """Add ``genes`` (overlapping) and ``promoter_of`` (TSS within pad) columns."""
    result = candidates.copy()
    if result.empty:
        result["genes"] = pd.Series(dtype=str)
        result["promoter_of"] = pd.Series(dtype=str)
        return result
    by_chrom = {chrom: group for chrom, group in genes.groupby("chrom")} if not genes.empty else {}
    overlapping, promoters = [], []
    for _, region in result.iterrows():
        chrom_genes = by_chrom.get(region["chrom"])
        if chrom_genes is None:
            overlapping.append("")
            promoters.append("")
            continue
        hit = chrom_genes[
            (chrom_genes["start"] < region["end"]) & (chrom_genes["end"] > region["start"])
        ]
        overlapping.append(",".join(sorted(hit["gene"].unique())[:limit]))
        prom = chrom_genes[
            (chrom_genes["tss"] - promoter_pad < region["end"])
            & (chrom_genes["tss"] + promoter_pad > region["start"])
        ]
        promoters.append(",".join(sorted(prom["gene"].unique())[:limit]))
    result["genes"] = overlapping
    result["promoter_of"] = promoters
    return result


def extract_to_regions(frame: pd.DataFrame) -> list[str]:
    """Convert a 0-based half-open BED frame into 1-based samtools region strings."""
    regions: list[str] = []
    for _, row in frame.iterrows():
        start = max(1, int(row["start"]) + 1)
        regions.append(f"{row['chrom']}:{start}-{int(row['end'])}")
    return regions


def write_bed3(frame: pd.DataFrame, path: str | Path) -> Path:
    """Write a plain 3-column BED (modkit --include-bed requires BED3 or BED6)."""
    output = Path(path)
    frame[["chrom", "start", "end"]].astype({"start": int, "end": int}).to_csv(
        output, sep="\t", header=False, index=False
    )
    return output


def panel_regions(
    gtf: str | Path, genes: list[str], promoter_pad: int, body_pad: int
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Build promoter + gene-body regions for a gene panel (mirrors script 06).

    Returns ``(named_regions, extract_intervals, missing_genes)`` where
    ``named_regions`` carries a ``name`` column formatted ``gene|promoter|body``
    for the targeted reader, and ``extract_intervals`` is a merged BED3 frame for
    restricting the pileup.
    """

    wanted = {symbol.strip().upper() for symbol in genes if symbol.strip()}
    found: dict[str, tuple[str, int, int, str]] = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = GENE_NAME.search(fields[8])
            if not match:
                continue
            key = GENE_ALIASES.get(match.group(1).upper(), match.group(1).upper())
            if key not in wanted:
                continue
            chrom, start, end, strand = fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]
            if key in found and (found[key][2] - found[key][1]) >= (end - start):
                continue
            found[key] = (chrom, start, end, strand)

    missing = sorted(wanted - set(found))
    named_rows: list[dict] = []
    extract_rows: list[tuple[str, int, int]] = []
    for gene, (chrom, start, end, strand) in sorted(found.items()):
        tss = start if strand == "+" else end
        named_rows.append({
            "chrom": chrom, "start": max(0, tss - promoter_pad), "end": tss + promoter_pad,
            "name": f"{gene}|promoter", "gene": gene, "region": "promoter",
        })
        named_rows.append({
            "chrom": chrom, "start": max(0, start - body_pad), "end": end + body_pad,
            "name": f"{gene}|body", "gene": gene, "region": "body",
        })
        extract_rows.append((chrom, max(0, start - body_pad - 1000), end + body_pad + 1000))

    named = pd.DataFrame(
        named_rows, columns=["chrom", "start", "end", "name", "gene", "region"]
    ).sort_values(["chrom", "start"]).reset_index(drop=True)

    extract_rows.sort()
    merged: list[list] = []
    for chrom, start, end in extract_rows:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    extract = pd.DataFrame(merged, columns=["chrom", "start", "end"])
    return named, extract, missing
