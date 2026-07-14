import gzip
import io
from pathlib import Path

import pytest

from core import references
from core.references import (
    ASSEMBLIES,
    Assembly,
    available_assemblies,
    is_prepared,
    prepare_assembly,
    prepared_paths,
    reference_cache_root,
)


class _FakeResponse:
    def __init__(self, payload: bytes, *, fail_after: int | None = None):
        self._stream = io.BytesIO(payload)
        self._fail_after = fail_after
        self._served = 0
        self.headers = {"Content-Length": str(len(payload))}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self, amt: int) -> bytes:
        if self._fail_after is not None and self._served >= self._fail_after:
            raise ConnectionError("simulated interruption")
        chunk = self._stream.read(amt)
        self._served += len(chunk)
        return chunk


def _gz(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb") as handle:
        handle.write(text.encode("utf-8"))
    return buffer.getvalue()


FASTA_GZ = _gz(">chr1\n" + "ACGT" * 20 + "\n")
GTF_GZ = _gz('chr1\tsource\tgene\t1\t100\t.\t+\t.\tgene_name "AAA";\n')
CPG_TXT_GZ = _gz("585\tchr1\t10\t50\tCpG:_5\t40\t35\n585\tchr1\t100\t160\tCpG:_9\n")


def _fixture_assembly(tmp_path: Path) -> Assembly:
    fasta = tmp_path / "hg.fa.gz"
    gtf = tmp_path / "hg.gtf.gz"
    cpg = tmp_path / "cpg.txt.gz"
    fasta.write_bytes(FASTA_GZ)
    gtf.write_bytes(GTF_GZ)
    cpg.write_bytes(CPG_TXT_GZ)
    return Assembly(
        key="hg38", label="Test build",
        fasta_url=fasta.as_uri(), gtf_url=gtf.as_uri(), cpg_url=cpg.as_uri(),
    )


def _payload_opener(mapping: dict[str, bytes], *, fail_after: int | None = None):
    calls: list[str] = []

    def opener(url: str):
        calls.append(url)
        return _FakeResponse(mapping[url], fail_after=fail_after)

    opener.calls = calls  # type: ignore[attr-defined]
    return opener


def test_manifests_cover_both_assemblies():
    assert set(available_assemblies()) == {"hg38", "hg19"}
    for key in available_assemblies():
        assembly = ASSEMBLIES[key]
        assert assembly.fasta_url.startswith("http")
        assert assembly.gtf_url.startswith("http")
        assert assembly.cpg_url.startswith("http")


def test_reference_cache_root_env(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("METHYL_TRIO_REFERENCE_CACHE", str(tmp_path / "cache"))
    assert reference_cache_root() == tmp_path / "cache"
    monkeypatch.delenv("METHYL_TRIO_REFERENCE_CACHE", raising=False)
    assert reference_cache_root() == Path.home() / ".cache" / "methyl-trio" / "references"


def test_prepare_assembly_downloads_prepares_and_reuses(monkeypatch, tmp_path: Path):
    assembly = _fixture_assembly(tmp_path)
    monkeypatch.setitem(ASSEMBLIES, "hg38", assembly)
    cache = tmp_path / "cache"

    mapping = {
        assembly.fasta_url: FASTA_GZ,
        assembly.gtf_url: GTF_GZ,
        assembly.cpg_url: CPG_TXT_GZ,
    }
    opener = _payload_opener(mapping)

    paths = prepare_assembly("hg38", opener=opener, cache_root=cache)

    assert paths["fasta"].read_text().startswith(">chr1")
    assert paths["fasta_index"].exists()
    assert paths["gtf"].exists()
    cpg_lines = paths["cpg_islands"].read_text().splitlines()
    assert cpg_lines[0].split("\t") == ["chr1", "10", "50", "CpG:_5"]
    assert is_prepared("hg38", cache_root=cache)
    assert len(opener.calls) == 3  # fasta, gtf, cpg each fetched once

    # A second preparation reuses the cache and performs no downloads.
    reuse_opener = _payload_opener(mapping)
    prepare_assembly("hg38", opener=reuse_opener, cache_root=cache)
    assert reuse_opener.calls == []


def test_interrupted_download_is_cleaned_up(monkeypatch, tmp_path: Path):
    assembly = _fixture_assembly(tmp_path)
    monkeypatch.setitem(ASSEMBLIES, "hg38", assembly)
    cache = tmp_path / "cache"

    mapping = {
        assembly.fasta_url: FASTA_GZ,
        assembly.gtf_url: GTF_GZ,
        assembly.cpg_url: CPG_TXT_GZ,
    }
    opener = _payload_opener(mapping, fail_after=0)

    with pytest.raises(ConnectionError):
        prepare_assembly("hg38", opener=opener, cache_root=cache)

    paths = prepared_paths("hg38", cache_root=cache)
    directory = paths["directory"]
    # No partial or finished artifacts should survive an interrupted download.
    if directory.exists():
        assert not any(directory.glob("*.part"))
        assert not paths["fasta"].exists()
    assert not is_prepared("hg38", cache_root=cache)
