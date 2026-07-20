"""High-level per-variant penetrance estimate combining the layer's parts."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

from penetrance.variant.betabinom import BetaBinomialPosterior, beta_binomial_posterior
from penetrance.variant.maxaf import penetrance_upper_bound_from_af


@dataclass(frozen=True)
class VariantPenetranceEstimate:
    """Result of the per-variant layer for a single variant."""

    point_estimate: float
    ci_low: float
    ci_high: float
    posterior: BetaBinomialPosterior
    af_upper_bound: Optional[float]
    confidence_tier: str
    data_weight: float
    prior: Tuple[float, float]
    notes: Dict[str, float] = field(default_factory=dict)


def _confidence_tier(n_observed: float, data_weight: float) -> str:
    """Confidence tier from the amount of direct data and prior dominance.

    ``data_weight`` is the share of the posterior's effective sample size that
    comes from observed carriers rather than the prior.
    """

    if n_observed >= 50 and data_weight >= 0.8:
        return "high"
    if n_observed >= 10 and data_weight >= 0.5:
        return "moderate"
    return "low"


def estimate_variant_penetrance(
    affected: float,
    unaffected: float,
    prior_alpha: float = 1.0,
    prior_beta: float = 1.0,
    allele_frequency: Optional[float] = None,
    prevalence: Optional[float] = None,
    f_allelic: float = 1.0,
    f_genetic: float = 1.0,
    level: float = 0.95,
    apply_af_bound: bool = True,
) -> VariantPenetranceEstimate:
    """Estimate penetrance for one variant.

    Combines the Beta-Binomial posterior (with a mechanism-aware prior) and,
    when allele frequency + prevalence are supplied, the Whiffin/Ware maximum
    credible allele-frequency upper bound. The bound only ever *lowers* the
    estimate - it can never manufacture penetrance that the counts do not
    support.
    """

    posterior = beta_binomial_posterior(affected, unaffected, prior_alpha, prior_beta)
    point = posterior.mean
    ci_low, ci_high = posterior.credible_interval(level)

    af_bound: Optional[float] = None
    notes: Dict[str, float] = {}
    if allele_frequency is not None and prevalence is not None:
        af_bound = penetrance_upper_bound_from_af(
            allele_frequency, prevalence, f_allelic, f_genetic
        )
        notes["af_upper_bound"] = af_bound
        if apply_af_bound:
            point = min(point, af_bound)
            ci_high = min(ci_high, af_bound)
            ci_low = min(ci_low, ci_high)

    n_obs = posterior.n_observed
    data_weight = n_obs / (n_obs + posterior.prior_strength) if n_obs + posterior.prior_strength > 0 else 0.0
    tier = _confidence_tier(n_obs, data_weight)

    return VariantPenetranceEstimate(
        point_estimate=float(point),
        ci_low=float(ci_low),
        ci_high=float(ci_high),
        posterior=posterior,
        af_upper_bound=af_bound,
        confidence_tier=tier,
        data_weight=float(data_weight),
        prior=(prior_alpha, prior_beta),
        notes=notes,
    )
