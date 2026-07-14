import gzip
from pathlib import Path

import pytest

from core import references


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    monkeypatch.setenv("METHYL_TRIO_REFERENCE_CACHE", str(tmp_path))
    return tmp_path


def test_assemblies_registered():
    assert set(references.ASSEMBLIES) == {"hg38", "hg19"}
    for spec in references.ASSEMBLIES.values():
        assert spec.fasta_url.startswith("https://")
        assert spec.gtf_url.endswith(".gtf.gz")


def test_bundle_paths_and_readiness(cache):
    bundle = references.bundle_paths("hg38")
    assert bundle.fasta.endswith("hg38/hg38.fa")
    assert not references.is_ready("hg38")
    directory = cache / "hg38"
    directory.mkdir()
    for path in (bundle.fasta, bundle.fasta_index, bundle.gtf, bundle.cpg_islands):
        Path(path).write_text("x")
    assert references.is_ready("hg38")


def test_ensure_unknown_assembly_raises():
    with pytest.raises(ValueError, match="Unknown assembly"):
        references.ensure_assembly("hg99")


def test_ensure_reuses_prepared_cache(cache, monkeypatch):
    bundle = references.bundle_paths("hg38")
    (cache / "hg38").mkdir()
    for path in (bundle.fasta, bundle.fasta_index, bundle.gtf, bundle.cpg_islands):
        Path(path).write_text("cached")

    def fail(*args, **kwargs):
        raise AssertionError("should not download when cache is ready")

    monkeypatch.setattr(references, "_download", fail)
    result = references.ensure_assembly("hg38")
    assert result == bundle


def test_interrupted_download_cleans_partial(cache, monkeypatch):
    destination = cache / "hg38" / "hg38.fa.gz"

    def broken(url, response=None):
        raise OSError("network dropped")

    class BrokenResponse:
        headers = {"Content-Length": "10"}

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self, size):
            raise OSError("network dropped")

    monkeypatch.setattr(references.urllib.request, "urlopen", lambda url: BrokenResponse())
    with pytest.raises(OSError):
        references._download("https://example/x", destination, lambda f, m: None, "x")
    assert not destination.with_suffix(destination.suffix + ".part").exists()
    assert not destination.exists()


def test_prepare_cpg_islands_extracts_bed_columns(cache):
    source = cache / "cpg.txt.gz"
    with gzip.open(source, "wt") as handle:
        handle.write("660\tchr1\t100\t200\tCpG:47\textra\n")
    destination = cache / "cpg.bed"
    references._prepare_cpg_islands(source, destination)
    assert destination.read_text().strip() == "chr1\t100\t200\tCpG:47"
