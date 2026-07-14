import pytest

from core.config import Affection, Relationship, Role, Sample, Sex, TrioConfig


def samples():
    return [
        Sample("P", "/p.bam", Sex.FEMALE, Role.PROBAND),
        Sample("M", "/m.bam", Sex.FEMALE, Role.RELATIVE),
        Sample("B", "/b.bam", Sex.MALE, Role.RELATIVE),
    ]


def test_comparisons_and_sex_chromosome_rules():
    config = TrioConfig(samples(), "/ref.fa")
    p_m, p_b, m_b = config.comparisons()
    assert p_m.valid_chromosome("chrX")
    assert not p_b.valid_chromosome("X")
    assert not m_b.valid_chromosome("chrY")
    assert not p_m.valid_chromosome("chrM")
    assert all(comparison.valid_chromosome("chr2") for comparison in config.comparisons())


def test_requires_exactly_one_proband_and_unique_labels():
    invalid = samples()
    invalid[1] = Sample("M", "/m.bam", Sex.FEMALE, Role.PROBAND)
    with pytest.raises(ValueError, match="Exactly one"):
        TrioConfig(invalid, "/ref.fa")

    duplicate = samples()
    duplicate[2] = Sample("M", "/b.bam", Sex.MALE, Role.RELATIVE)
    with pytest.raises(ValueError, match="unique"):
        TrioConfig(duplicate, "/ref.fa")


def test_optional_metadata_drives_evidence_status():
    trio = samples()
    trio[0] = Sample(
        "P", "/p.bam", Sex.FEMALE, Role.PROBAND,
        affection=Affection.AFFECTED, tissue="blood", batch="run-1",
    )
    trio[1] = Sample(
        "M", "/m.bam", Sex.FEMALE, Role.RELATIVE,
        Relationship.MOTHER, Affection.AFFECTED, "blood", "run-1",
    )
    trio[2] = Sample(
        "F", "/f.bam", Sex.MALE, Role.RELATIVE,
        Relationship.FATHER, Affection.UNAFFECTED, "blood", "run-1",
    )
    config = TrioConfig(trio, "/ref.fa", phased_vcf="/family.vcf.gz")
    assert config.analysis_design() == "phenotype_segregation"
    assert config.evidence_status() == {
        "phenotype": "phenotype_segregation",
        "parent_of_origin": "inputs_available",
        "mqtl": "phased_vcf_available",
        "tissue": "matched",
        "batch": "matched",
    }
