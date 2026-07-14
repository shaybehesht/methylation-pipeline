"""One-time cached download and preparation of managed reference assemblies."""
from __future__ import annotations

import gzip
import os
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

Progress = Callable[[float, str], None]


@dataclass(frozen=True)
class Assembly:
    key: str
    label: str
    fasta_url: str
    gtf_url: str
    cpg_islands_url: str
    ucsc_db: str


ASSEMBLIES: dict[str, Assembly] = {
    "hg38": Assembly(
        key="hg38",
        label="GRCh38 / hg38",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
        cpg_islands_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz",
        ucsc_db="hg38",
    ),
    "hg19": Assembly(
        key="hg19",
        label="GRCh37 / hg19",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh37_mapping/gencode.v49lift37.annotation.gtf.gz",
        cpg_islands_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz",
        ucsc_db="hg19",
    ),
}


@dataclass(frozen=True)
class ReferenceBundle:
    assembly: str
    fasta: str
    fasta_index: str
    gtf: str
    cpg_islands: str


def cache_root() -> Path:
    configured = os.environ.get("METHYL_TRIO_REFERENCE_CACHE")
    root = Path(configured) if configured else Path.home() / ".cache" / "methyl-trio" / "references"
    return root.expanduser().resolve()


def _download(url: str, destination: Path, notify: Progress, message: str) -> None:
    """Stream a URL to ``destination`` atomically via a temporary file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with urllib.request.urlopen(url) as response:  # noqa: S310 - fixed manifest URLs
            total = int(response.headers.get("Content-Length", 0))
            read = 0
            with temporary.open("wb") as handle:
                while True:
                    chunk = response.read(1 << 20)
                    if not chunk:
                        break
                    handle.write(chunk)
                    read += len(chunk)
                    fraction = read / total if total else 0.0
                    notify(min(fraction, 0.99), f"{message} ({read // (1 << 20)} MB)")
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _decompress(source: Path, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with gzip.open(source, "rb") as compressed, temporary.open("wb") as handle:
            shutil.copyfileobj(compressed, handle)
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _prepare_cpg_islands(source_gz: Path, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with gzip.open(source_gz, "rt") as handle, temporary.open("w", encoding="utf-8") as out:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                out.write(f"{fields[1]}\t{fields[2]}\t{fields[3]}\t{fields[4]}\n")
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def bundle_paths(assembly: str) -> ReferenceBundle:
    directory = cache_root() / assembly
    return ReferenceBundle(
        assembly=assembly,
        fasta=str(directory / f"{assembly}.fa"),
        fasta_index=str(directory / f"{assembly}.fa.fai"),
        gtf=str(directory / "annotation.gtf.gz"),
        cpg_islands=str(directory / "cpg_islands.bed"),
    )


def is_ready(assembly: str) -> bool:
    bundle = bundle_paths(assembly)
    return all(
        Path(path).exists()
        for path in (bundle.fasta, bundle.fasta_index, bundle.gtf, bundle.cpg_islands)
    )


def ensure_assembly(assembly: str, progress: Progress | None = None) -> ReferenceBundle:
    """Download, decompress, index, and cache an assembly. Idempotent and offline once cached."""
    if assembly not in ASSEMBLIES:
        raise ValueError(f"Unknown assembly: {assembly}")
    notify = progress or (lambda fraction, message: None)
    spec = ASSEMBLIES[assembly]
    bundle = bundle_paths(assembly)
    directory = cache_root() / assembly
    directory.mkdir(parents=True, exist_ok=True)

    if is_ready(assembly):
        notify(1.0, "Reference already prepared")
        return bundle

    fasta = Path(bundle.fasta)
    if not fasta.exists():
        fasta_gz = directory / f"{assembly}.fa.gz"
        if not fasta_gz.exists():
            notify(0.05, "Downloading reference FASTA")
            _download(spec.fasta_url, fasta_gz, notify, "Downloading reference FASTA")
        notify(0.5, "Decompressing reference FASTA")
        _decompress(fasta_gz, fasta)
        fasta_gz.unlink(missing_ok=True)

    if not Path(bundle.fasta_index).exists():
        notify(0.7, "Indexing reference FASTA")
        import pysam

        pysam.faidx(str(fasta))

    if not Path(bundle.gtf).exists():
        notify(0.85, "Downloading gene annotation")
        _download(spec.gtf_url, Path(bundle.gtf), notify, "Downloading gene annotation")

    if not Path(bundle.cpg_islands).exists():
        notify(0.95, "Downloading CpG islands")
        cpg_gz = directory / "cpgIslandExt.txt.gz"
        _download(spec.cpg_islands_url, cpg_gz, notify, "Downloading CpG islands")
        _prepare_cpg_islands(cpg_gz, Path(bundle.cpg_islands))
        cpg_gz.unlink(missing_ok=True)

    notify(1.0, "Reference ready")
    return bundle
