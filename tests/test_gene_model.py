import numpy as np
import pandas as pd
import pytest

from penetrance.features.matrix import FEATURE_COLUMNS, build_feature_matrix
from penetrance.features.mechanism import infer_dn_gof_score, infer_mechanism_class
from penetrance.gene_model.cv import GeneFamilyKFold, family_aware_split
from penetrance.gene_model.prior import (
    beta_prior_from_mean_var,
    propensity_to_beta_prior,
)
from penetrance.gene_model.model import train_gene_model
from penetrance.labels.loader import load_gene_labels


@pytest.fixture(scope="module")
def genes():
    return load_gene_labels()


def test_mechanism_inference_tubulin_vs_brca(genes):
    scores = infer_dn_gof_score(genes.set_index("gene"))
    assert scores["TUBA1A"] > 0.7  # dominant-negative / GoF
    assert scores["BRCA1"] < 0.3  # LoF / haploinsufficiency
    assert scores["MYH7"] > scores["MYBPC3"]  # DN vs HI within sarcomere


def test_mechanism_class_labels(genes):
    classes = infer_mechanism_class(genes.set_index("gene"))
    assert classes["TUBA1A"] == "DN_GoF"
    assert classes["BRCA1"] == "LoF_HI"


def test_feature_matrix_shape_and_columns(genes):
    X, meta = build_feature_matrix(genes)
    assert list(X.columns) == FEATURE_COLUMNS
    assert len(X) == len(genes)
    assert not X.isna().any().any()
    assert set(meta.columns) == {"gene", "gene_family"}


def test_family_aware_cv_has_no_family_leakage(genes):
    groups = genes["gene_family"].to_numpy()
    splitter = GeneFamilyKFold(n_splits=5, random_state=1)
    seen_test = 0
    for train_idx, test_idx in splitter.split(genes, genes["penetrance"], groups):
        train_fams = set(groups[train_idx])
        test_fams = set(groups[test_idx])
        assert train_fams.isdisjoint(test_fams)
        seen_test += len(test_idx)
    assert seen_test == len(genes)  # every gene tested exactly once


def test_family_aware_split_holds_out_named_family(genes):
    groups = genes["gene_family"].to_numpy()
    train_idx, test_idx = family_aware_split(groups, ["tubulin"])
    assert set(groups[test_idx]) == {"tubulin"}
    assert "tubulin" not in set(groups[train_idx])


def test_prior_mapping_mean_and_strength():
    prior = propensity_to_beta_prior(0.7, strength=10)
    assert prior.mean == pytest.approx(0.7, abs=1e-3)
    assert prior.strength == pytest.approx(10.0)


def test_prior_from_mean_var_matches_moments():
    m, v = 0.4, 0.02
    prior = beta_prior_from_mean_var(m, v)
    assert prior.mean == pytest.approx(m, abs=1e-3)
    assert prior.variance == pytest.approx(v, rel=0.05)


def test_gene_model_trains_and_separates_mechanism(genes):
    model, cv = train_gene_model(genes, n_splits=5, random_state=0)
    # Out-of-fold rank correlation should be clearly positive.
    assert cv.metrics["spearman"] > 0.4
    X, meta = build_feature_matrix(genes)
    preds = pd.Series(model.predict(X), index=meta["gene"].values)
    tubulins = preds[genes.set_index("gene").loc[preds.index, "gene_family"].values == "tubulin"]
    hboc = preds[genes.set_index("gene").loc[preds.index, "gene_family"].values == "hboc"]
    assert tubulins.mean() > hboc.mean()


def test_gene_model_beta_prior_fixed_strength(genes):
    model, _ = train_gene_model(genes, prior_strength=12.0)
    X, _ = build_feature_matrix(genes.head(3))
    priors = model.beta_prior_for(X)
    for p in priors:
        assert p.strength == pytest.approx(12.0)
        assert 0 <= p.mean <= 1
