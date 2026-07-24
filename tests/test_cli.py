import pytest

from core import cli
from core.config import Affection, Relationship, Role, Sex


BASE_ARGS = [
    "--proband-bam", "proband.bam", "--proband-sex", "F", "--proband-affection", "affected",
    "--relative1-bam", "mother.bam", "--relative1-sex", "F",
    "--relative1-relationship", "mother", "--relative1-affection", "unaffected",
    "--relative2-bam", "father.bam", "--relative2-sex", "M",
    "--relative2-relationship", "father",
    "--reference-fasta", "hg38.fa",
    "--output-dir", "out",
]


def _config(extra):
    args = cli.build_parser().parse_args(BASE_ARGS + extra)
    return cli.build_config(args)


def test_build_config_targeted():
    config = _config(["--gtf", "g.gtf", "--mode", "targeted", "--genes", "MECP2", "UBE3A"])
    assert len(config.samples) == 3
    proband, r1, r2 = config.samples
    assert proband.role == Role.PROBAND and proband.sex == Sex.FEMALE
    assert proband.affection == Affection.AFFECTED
    assert r1.role == Role.RELATIVE and r1.relationship == Relationship.MOTHER
    assert r1.affection == Affection.UNAFFECTED
    assert r2.sex == Sex.MALE and r2.relationship == Relationship.FATHER
    assert r2.affection is None
    assert config.regions.mode == "targeted"
    assert config.regions.genes == ["MECP2", "UBE3A"]
    assert config.reference_fasta == "hg38.fa"
    assert config.output_dir == "out"
    assert config.combine_strands is True
    assert config.modified_bases == ["5mC"]


def test_build_config_chromosomes_and_options():
    config = _config([
        "--mode", "chromosomes", "--chromosomes", "chr1", "chr11",
        "--no-combine-strands", "--modified-base", "5mC", "--modified-base", "5hmC",
    ])
    assert config.regions.mode == "chromosomes"
    assert config.regions.chromosomes == ["chr1", "chr11"]
    assert config.combine_strands is False
    assert config.modified_bases == ["5mC", "5hmC"]
    assert config.platform == "auto"


def test_build_config_platform_ont_and_pacbio():
    ont = _config(["--mode", "whole_genome", "--platform", "ont"])
    assert ont.platform == "ont"
    assert ont.combine_strands is True
    pacbio = _config(["--mode", "whole_genome", "--platform", "pacbio"])
    assert pacbio.platform == "pacbio"


def test_build_config_whole_genome_defaults():
    config = _config(["--mode", "whole_genome"])
    assert config.regions.mode == "whole_genome"


def test_default_label_used_when_absent():
    config = _config(["--mode", "whole_genome"])
    assert [s.label for s in config.samples] == ["proband", "relative1", "relative2"]


def test_custom_labels():
    config = _config([
        "--mode", "whole_genome",
        "--proband-label", "BH1", "--relative1-label", "Mom", "--relative2-label", "Dad",
    ])
    assert [s.label for s in config.samples] == ["BH1", "Mom", "Dad"]


def test_targeted_without_genes_is_rejected():
    args = cli.build_parser().parse_args(BASE_ARGS + ["--mode", "targeted"])
    # RegionConfig enforces "targeted needs genes" (build_config surfaces ValueError;
    # main() turns this into a parser error / non-zero exit).
    with pytest.raises(ValueError, match="target gene"):
        cli.build_config(args)


def test_main_rejects_targeted_without_genes():
    with pytest.raises(SystemExit):
        cli.main(BASE_ARGS + ["--mode", "targeted"])


def test_parse_thresholds_valid_and_types():
    result = cli.parse_thresholds(["min_sites=8", "alpha=0.01", "filter_threshold=0.75"])
    assert result["min_sites"] == 8 and isinstance(result["min_sites"], int)
    assert result["alpha"] == 0.01
    assert result["filter_threshold"] == 0.75
    # unspecified keys keep their defaults
    assert "null_percentile" in result


def test_parse_thresholds_rejects_unknown_and_malformed():
    with pytest.raises(ValueError, match="Unknown threshold"):
        cli.parse_thresholds(["nope=1"])
    with pytest.raises(ValueError, match="KEY=VALUE"):
        cli.parse_thresholds(["min_sites"])


def test_set_threshold_flows_into_config():
    config = _config(["--mode", "whole_genome", "--set-threshold", "min_sites=7"])
    assert config.thresholds["min_sites"] == 7


# --- reference resolution (auto-download) ---

NO_REF_ARGS = [
    "--proband-bam", "proband.bam", "--proband-sex", "F",
    "--relative1-bam", "mother.bam", "--relative1-sex", "F",
    "--relative2-bam", "father.bam", "--relative2-sex", "M",
    "--output-dir", "out", "--mode", "whole_genome",
]


def test_reference_fasta_is_optional_in_parser():
    args = cli.build_parser().parse_args(NO_REF_ARGS)
    assert args.reference_fasta == ""


def test_resolve_reference_noop_when_fasta_provided():
    args = cli.build_parser().parse_args(BASE_ARGS + ["--mode", "whole_genome"])
    cli.resolve_reference(args)  # must not touch the network
    assert args.reference_fasta == "hg38.fa"


def test_resolve_reference_requires_fasta_or_assembly():
    args = cli.build_parser().parse_args(NO_REF_ARGS)
    with pytest.raises(ValueError, match="assembly"):
        cli.resolve_reference(args)


def test_resolve_reference_downloads_from_assembly(monkeypatch, tmp_path):
    import core.references as refs

    def fake_prepare(key, progress=None):
        assert key == "hg38"
        return {
            "fasta": tmp_path / "genome.fa",
            "gtf": tmp_path / "gencode.gtf.gz",
            "cpg_islands": tmp_path / "cpg.bed",
        }

    monkeypatch.setattr(refs, "prepare_assembly", fake_prepare)
    args = cli.build_parser().parse_args(NO_REF_ARGS + ["--assembly", "hg38"])
    cli.resolve_reference(args)
    assert args.reference_fasta == str(tmp_path / "genome.fa")
    assert args.gtf == str(tmp_path / "gencode.gtf.gz")
    assert args.cpg_islands == str(tmp_path / "cpg.bed")
    # a resolved config can now be built
    config = cli.build_config(args)
    assert config.reference_fasta == str(tmp_path / "genome.fa")
