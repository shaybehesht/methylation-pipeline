import pytest

from penetrance.adapters.literature import LiteratureCountAdapter, LiteratureRecord
from penetrance.adapters.clinvar_gnomad import FrequencyCountAdapter
from penetrance.pipeline import PenetrancePipeline
from penetrance.validation import (
    gene_model_calibration,
    mechanism_separation,
    sparse_prior_benefit,
)


@pytest.fixture(scope="module")
def pipeline():
    return PenetrancePipeline().fit()


def test_pipeline_gene_prior_orders_mechanism(pipeline):
    tub = pipeline.gene_prior("TUBA1A")
    brca = pipeline.gene_prior("BRCA1")
    assert tub.mean > brca.mean
    assert pipeline.gene_prior("BRCA1", flat=True).as_tuple() == (1.0, 1.0)


def test_variant_estimate_flat_vs_mechanism(pipeline):
    est = pipeline.estimate_variant("SDHB:p.Ile127Ser", use_gene_prior=True)
    flat = pipeline.estimate_variant("SDHB:p.Ile127Ser", use_gene_prior=False)
    assert est is not None and flat is not None
    assert 0 <= est.point_estimate <= 1
    # The mechanism prior should tighten the credible interval vs the flat prior.
    mech_width = est.ci_high - est.ci_low
    flat_width = flat.ci_high - flat.ci_low
    assert mech_width <= flat_width


def test_af_bound_caps_hemochromatosis(pipeline):
    est = pipeline.estimate_variant("HFE:p.Cys282Tyr")
    assert est.point_estimate < 0.2  # classic incomplete penetrance


def test_mechanism_separation_positive(pipeline):
    sep = mechanism_separation(pipeline)
    assert sep["separation"] > 0.2


def test_sparse_prior_beats_flat(pipeline):
    res = sparse_prior_benefit(pipeline, subsample_sizes=(1, 3), n_trials=150)
    improvement = res.improvement()
    # Mechanism prior must reduce error for the sparsest case.
    assert improvement.loc[res.summary["k"] == 1].iloc[0] > 0
    assert (res.summary["mechanism_mae"] <= res.summary["flat_mae"]).all()


def test_gene_model_calibration_runs(pipeline):
    calib = gene_model_calibration(pipeline)
    assert 0 <= calib["expected_calibration_error"] <= 1
    assert calib["metrics"]["spearman"] > 0.3


def test_pipeline_with_multiple_adapters():
    records = [
        LiteratureRecord(
            doc_id="PMID:99",
            gene="SDHB",
            text="For p.Ile127Ser, 5 of 8 carriers were affected in the extended pedigree.",
        )
    ]
    pipe = PenetrancePipeline(
        adapters=[FrequencyCountAdapter(), LiteratureCountAdapter(records)]
    ).fit()
    counts = pipe.carrier_counts("SDHB:p.Ile127Ser")
    # Frequency (9 affected / 17 unaffected) + literature (5 affected / 3 unaffected).
    assert counts.affected == 14
    assert counts.unaffected == 20
    assert "literature" in counts.source
