import math

import pytest

from penetrance.variant.betabinom import beta_binomial_posterior
from penetrance.variant.casecontrol import (
    case_control_penetrance,
    case_control_penetrance_from_counts,
)
from penetrance.variant.maxaf import max_credible_af, penetrance_upper_bound_from_af
from penetrance.variant.estimator import estimate_variant_penetrance


def test_beta_binomial_posterior_mean():
    post = beta_binomial_posterior(affected=8, unaffected=2, prior_alpha=1, prior_beta=1)
    # (1+8)/(1+1+8+2) = 9/12
    assert post.mean == pytest.approx(9 / 12)
    assert post.alpha == 9 and post.beta == 3


def test_beta_binomial_flat_prior_matches_empirical_in_limit():
    post = beta_binomial_posterior(900, 100, 1, 1)
    assert post.mean == pytest.approx(0.9, abs=1e-3)


def test_credible_interval_contains_mean_and_narrows_with_data():
    small = beta_binomial_posterior(3, 1, 1, 1)
    large = beta_binomial_posterior(300, 100, 1, 1)
    lo_s, hi_s = small.credible_interval()
    lo_l, hi_l = large.credible_interval()
    assert lo_s <= small.mean <= hi_s
    assert (hi_l - lo_l) < (hi_s - lo_s)


def test_prior_pulls_sparse_estimate():
    # One affected carrier, strong prior toward low penetrance.
    flat = beta_binomial_posterior(1, 0, 1, 1).mean
    low_prior = beta_binomial_posterior(1, 0, 2, 10).mean
    assert low_prior < flat


def test_invalid_inputs():
    with pytest.raises(ValueError):
        beta_binomial_posterior(-1, 0)
    with pytest.raises(ValueError):
        beta_binomial_posterior(1, 1, prior_alpha=0)


def test_case_control_formula_known_value():
    # freq cases 0.5, controls 0.0, any prevalence -> penetrance 1.0
    assert case_control_penetrance(0.5, 0.0, 0.01) == pytest.approx(1.0)
    # equal frequencies -> penetrance == prevalence
    assert case_control_penetrance(0.2, 0.2, 0.05) == pytest.approx(0.05)


def test_case_control_from_counts_ci_brackets_point():
    est = case_control_penetrance_from_counts(30, 100, 2, 1000, prevalence=0.01, n_boot=5000)
    assert est.ci_low <= est.penetrance <= est.ci_high
    assert 0 <= est.ci_low <= est.ci_high <= 1


def test_max_credible_af_and_inverse_are_consistent():
    prev, pen = 0.001, 0.5
    af = max_credible_af(prev, pen, f_allelic=0.5, f_genetic=0.5)
    # Inverting at that AF should recover the penetrance.
    bound = penetrance_upper_bound_from_af(af, prev, f_allelic=0.5, f_genetic=0.5)
    assert bound == pytest.approx(pen, rel=1e-6)


def test_af_bound_is_clamped_and_monotone():
    high = penetrance_upper_bound_from_af(1e-6, 0.001)
    low = penetrance_upper_bound_from_af(0.05, 0.001)
    assert high == 1.0  # rare allele -> uninformative bound
    assert low < high  # common allele -> penetrance capped below 1


def test_estimator_applies_af_bound_only_downward():
    # A common allele in a rare disease must cap penetrance low.
    est = estimate_variant_penetrance(
        affected=50, unaffected=5, allele_frequency=0.05, prevalence=0.001
    )
    assert est.af_upper_bound is not None
    assert est.point_estimate <= est.af_upper_bound + 1e-9
    assert est.point_estimate < 0.5  # counts alone would say ~0.9


def test_estimator_confidence_tiers():
    high = estimate_variant_penetrance(90, 10, prior_alpha=1, prior_beta=1)
    low = estimate_variant_penetrance(1, 0, prior_alpha=6, prior_beta=6)
    assert high.confidence_tier == "high"
    assert low.confidence_tier == "low"
