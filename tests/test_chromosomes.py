from pathlib import Path

import pytest

from core.chromosomes import dmr_reference_contig, pileup_contig_name, resolve_reference_contig


def test_pileup_contig_name_reads_first_data_row(tmp_path: Path):
    path = tmp_path / "pileup.bed"
    path.write_text("#comment\n14\t100\t101\tm\n", encoding="utf-8")
    assert pileup_contig_name(path) == "14"


def test_resolve_reference_contig_prefers_chr_prefix(tmp_path: Path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr14\nACGT\n", encoding="utf-8")
    assert resolve_reference_contig(str(fasta), "chr14") == "chr14"
    assert resolve_reference_contig(str(fasta), "14") == "chr14"


def test_dmr_reference_contig_renames_fasta_header_to_match_pileup(tmp_path: Path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr14\nACGT\n", encoding="utf-8")
    pileup = tmp_path / "sample.bed"
    pileup.write_text("14\t1\t2\tm\n", encoding="utf-8")
    ref_contig, header = dmr_reference_contig(str(fasta), "chr14", [str(pileup)])
    assert ref_contig == "chr14"
    assert header == "14"


def test_dmr_reference_contig_rejects_mismatched_pileups(tmp_path: Path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr14\nACGT\n", encoding="utf-8")
    left = tmp_path / "left.bed"
    right = tmp_path / "right.bed"
    left.write_text("14\t1\t2\tm\n", encoding="utf-8")
    right.write_text("chr14\t1\t2\tm\n", encoding="utf-8")
    with pytest.raises(ValueError, match="disagree"):
        dmr_reference_contig(str(fasta), "chr14", [str(left), str(right)])
