import pytest

from penetrance.adapters.base import combine_counts, CarrierCounts
from penetrance.adapters.clinvar_gnomad import FrequencyCountAdapter
from penetrance.adapters.literature import (
    LiteratureCountAdapter,
    LiteratureRecord,
    normalize_variant_tokens,
    regex_carrier_extractor,
)


def test_frequency_adapter_fetches_counts():
    adapter = FrequencyCountAdapter()
    counts = adapter.fetch("HFE:p.Cys282Tyr")
    assert counts is not None
    assert counts.affected == 55
    assert counts.unaffected == 480
    assert counts.allele_frequency == pytest.approx(0.06)
    assert adapter.fetch("does:not:exist") is None


def test_variant_token_normalization_protein_forms():
    tokens = normalize_variant_tokens(
        hgvs="NM_006009.4:c.1205G>A", variant_id="TUBA1A:p.Arg402His"
    )
    assert "r402h" in tokens  # 1-letter form
    assert "p.r402h" in tokens
    assert any(t.startswith("c.1205") for t in tokens)


def test_literature_extractor_and_variant_gate():
    records = [
        LiteratureRecord(
            doc_id="PMID:1",
            gene="TUBA1A",
            text=(
                "We studied the p.Arg402His variant. "
                "In this cohort, 12 of 15 carriers were affected. "
                "An unrelated variant p.Gly99Ser showed 1 of 40 carriers affected."
            ),
            metadata={"phenotype": "lissencephaly"},
        )
    ]
    adapter = LiteratureCountAdapter(records)
    counts = adapter.fetch("TUBA1A:p.Arg402His")
    assert counts is not None
    # Only the sentence mentioning the target variant should be counted.
    assert counts.affected == 12
    assert counts.unaffected == 3
    assert counts.provenance and "PMID:1" in counts.provenance[0]


def test_literature_gene_gate_rejects_other_gene():
    records = [
        LiteratureRecord(
            doc_id="PMID:2",
            gene="BRCA1",
            text="p.Arg402His had 9 of 10 carriers affected.",
        )
    ]
    # Target variant belongs to TUBA1A; BRCA1 document must be skipped.
    adapter = LiteratureCountAdapter(records)
    assert adapter.fetch("TUBA1A:p.Arg402His") is None


def test_regex_extractor_affected_unaffected_pattern():
    rec = LiteratureRecord(
        doc_id="D",
        gene="X",
        text="For c.1205G>A there were 4 affected and 6 unaffected carriers.",
    )
    results = regex_carrier_extractor(rec, ["c.1205g>a"])
    assert len(results) == 1
    assert results[0].affected == 4 and results[0].unaffected == 6


def test_combine_counts_sums_and_keeps_af():
    a = CarrierCounts("v", 3, 5, "clinvar_gnomad", allele_frequency=1e-4, provenance=["p1"])
    b = CarrierCounts("v", 4, 1, "literature", provenance=["p2"])
    merged = combine_counts([a, b])
    assert merged.affected == 7 and merged.unaffected == 6
    assert merged.allele_frequency == pytest.approx(1e-4)
    assert merged.provenance == ["p1", "p2"]
    assert "clinvar_gnomad" in merged.source and "literature" in merged.source
