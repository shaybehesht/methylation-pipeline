"""Managed hg38/hg19 reference assemblies with a persistent local cache.

The first time an assembly is requested its FASTA, GENCODE annotation, and UCSC
CpG-island track are streamed from public mirrors, written atomically, and
prepared (decompressed and indexed). Every later run reuses the cache and needs
no network access. The cache location is controlled by
``METHYL_TRIO_REFERENCE_CACHE`` so it can be persisted through a host-mounted
directory or Docker volume instead of being baked into the image.
"""
from __future__ import annotations

import gzip
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, ContextManager, Protocol
from urllib.request import urlopen

import pysam

Progress = Callable[[float, str], None]
_CHUNK = 1 << 20  # 1 MiB streaming chunks keep memory flat on large FASTAs.


class _Response(Protocol):
    headers: object

    def read(self, amt: int) -> bytes: ...


Opener = Callable[[str], ContextManager[_Response]]


@dataclass(frozen=True)
class Assembly:
    """Download manifest for a single genome build."""

    key: str
    label: str
    fasta_url: str
    gtf_url: str
    cpg_url: str


ASSEMBLIES: dict[str, Assembly] = {
    "hg38": Assembly(
        key="hg38",
        label="GRCh38 / hg38",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        gtf_url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_49/gencode.v49.annotation.gtf.gz"
        ),
        cpg_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz",
    ),
    "hg19": Assembly(
        key="hg19",
        label="GRCh37 / hg19",
        fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        gtf_url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_49/GRCh37_mapping/gencode.v49lift37.annotation.gtf.gz"
        ),
        cpg_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz",
    ),
}


def available_assemblies() -> list[str]:
    return list(ASSEMBLIES)


def get_assembly(key: str) -> Assembly:
    try:
        return ASSEMBLIES[key]
    except KeyError:
        raise ValueError(f"Unknown assembly: {key!r}. Choose one of {available_assemblies()}.")


def reference_cache_root() -> Path:
    """Return the root directory for cached assemblies."""

    configured = os.environ.get("METHYL_TRIO_REFERENCE_CACHE")
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".cache" / "methyl-trio" / "references"


def assembly_dir(key: str, cache_root: str | Path | None = None) -> Path:
    root = Path(cache_root) if cache_root is not None else reference_cache_root()
    return root / key


def prepared_paths(key: str, cache_root: str | Path | None = None) -> dict[str, Path]:
    """Return the canonical prepared artifact paths for an assembly."""

    directory = assembly_dir(key, cache_root)
    return {
        "assembly": Path(key),
        "directory": directory,
        "fasta": directory / "genome.fa",
        "fasta_index": directory / "genome.fa.fai",
        "gtf": directory / "gencode.annotation.gtf.gz",
        "cpg_islands": directory / "cpg_islands.bed",
    }


def is_prepared(key: str, cache_root: str | Path | None = None) -> bool:
    """Return whether every prepared artifact already exists in the cache."""

    paths = prepared_paths(key, cache_root)
    return all(
        paths[name].exists()
        for name in ("fasta", "fasta_index", "gtf", "cpg_islands")
    )


def _stream_download(
    url: str, dest: Path, *, opener: Opener, progress: Progress | None, label: str
) -> Path:
    """Stream ``url`` to ``dest`` atomically, cleaning up on any interruption."""

    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    try:
        with opener(url) as response:
            header_value = response.headers.get("Content-Length") if response.headers else None
            total = int(header_value) if header_value else 0
            downloaded = 0
            with partial.open("wb") as handle:
                while True:
                    chunk = response.read(_CHUNK)
                    if not chunk:
                        break
                    handle.write(chunk)
                    downloaded += len(chunk)
                    if progress:
                        fraction = downloaded / total if total else 0.0
                        megabytes = downloaded / 1_000_000
                        progress(fraction, f"Downloading {label} ({megabytes:.0f} MB)")
        os.replace(partial, dest)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    return dest


def _decompress(src: Path, dest: Path) -> Path:
    """Gunzip ``src`` into ``dest`` atomically, cleaning up on failure."""

    partial = dest.with_name(dest.name + ".part")
    try:
        with gzip.open(src, "rb") as compressed, partial.open("wb") as out:
            shutil.copyfileobj(compressed, out, _CHUNK)
        os.replace(partial, dest)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    return dest


def index_fasta(fasta: str | Path) -> Path:
    """Build a ``.fai`` index for ``fasta`` using pysam."""

    pysam.faidx(str(fasta))
    return Path(str(fasta) + ".fai")


def prepare_cpg_islands(source_txt_gz: Path, dest_bed: Path) -> Path:
    """Convert a UCSC ``cpgIslandExt.txt.gz`` file into a 4-column BED.

    UCSC columns are ``bin, chrom, chromStart, chromEnd, name, ...``; the
    prepared BED keeps ``chrom, start, end, name`` to match the pipeline scope
    reader.
    """

    partial = dest_bed.with_name(dest_bed.name + ".part")
    try:
        with gzip.open(source_txt_gz, "rt") as handle, partial.open("w", encoding="utf-8") as out:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                out.write("\t".join((fields[1], fields[2], fields[3], fields[4])) + "\n")
        os.replace(partial, dest_bed)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    return dest_bed


def prepare_assembly(
    key: str,
    *,
    progress: Progress | None = None,
    opener: Opener = urlopen,
    cache_root: str | Path | None = None,
) -> dict[str, Path]:
    """Download and prepare ``key`` if needed, returning prepared artifact paths.

    Each artifact is skipped when already cached, so repeated calls are cheap
    and fully offline once the first preparation succeeds.
    """

    assembly = get_assembly(key)
    directory = assembly_dir(key, cache_root)
    directory.mkdir(parents=True, exist_ok=True)
    paths = prepared_paths(key, cache_root)

    def notify(fraction: float, message: str) -> None:
        if progress:
            progress(fraction, message)

    if not paths["fasta"].exists():
        notify(0.0, f"Fetching {assembly.label} FASTA")
        compressed = directory / "genome.fa.gz"
        _stream_download(
            assembly.fasta_url, compressed, opener=opener, progress=progress,
            label=f"{assembly.label} FASTA",
        )
        notify(0.45, "Decompressing FASTA")
        _decompress(compressed, paths["fasta"])
        compressed.unlink(missing_ok=True)

    if not paths["fasta_index"].exists():
        notify(0.6, "Indexing FASTA with pysam")
        index_fasta(paths["fasta"])

    if not paths["gtf"].exists():
        notify(0.7, f"Fetching {assembly.label} GENCODE annotation")
        _stream_download(
            assembly.gtf_url, paths["gtf"], opener=opener, progress=progress,
            label=f"{assembly.label} GTF",
        )

    if not paths["cpg_islands"].exists():
        notify(0.9, "Fetching CpG islands")
        compressed = directory / "cpgIslandExt.txt.gz"
        _stream_download(
            assembly.cpg_url, compressed, opener=opener, progress=progress,
            label=f"{assembly.label} CpG islands",
        )
        prepare_cpg_islands(compressed, paths["cpg_islands"])
        compressed.unlink(missing_ok=True)

    notify(1.0, f"{assembly.label} ready")
    return paths
