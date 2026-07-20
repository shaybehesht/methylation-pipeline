"""Map a gene-level penetrance propensity to a Beta prior ``(alpha, beta)``.

The gene model outputs a continuous propensity ``m in [0, 1]`` (the predicted
penetrance tendency). We turn it into the prior for the per-variant
Beta-Binomial layer. Two mappings are provided:

* :func:`propensity_to_beta_prior` - fix the mean to ``m`` and the concentration
  (prior "strength" / pseudo-count) to a chosen value ``kappa``. Larger ``kappa``
  means we trust the gene prediction more and it takes more carrier
  observations to overturn it.
* :func:`beta_prior_from_mean_var` - method-of-moments prior from a predicted
  mean and variance (use the model's per-gene predictive uncertainty).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

_EPS = 1e-6


@dataclass(frozen=True)
class BetaPrior:
    alpha: float
    beta: float

    @property
    def mean(self) -> float:
        return self.alpha / (self.alpha + self.beta)

    @property
    def strength(self) -> float:
        """Concentration = effective prior sample size (pseudo-count)."""

        return self.alpha + self.beta

    @property
    def variance(self) -> float:
        a, b = self.alpha, self.beta
        return (a * b) / ((a + b) ** 2 * (a + b + 1))

    def as_tuple(self) -> Tuple[float, float]:
        return (self.alpha, self.beta)


def propensity_to_beta_prior(propensity: float, strength: float = 10.0) -> BetaPrior:
    """Beta prior with mean ``propensity`` and concentration ``strength``.

    ``strength`` is the pseudo-count: with ``strength=10`` the gene prior carries
    the weight of ~10 observed carriers, so a variant with many more real
    carriers will dominate it, while a singleton variant leans on the prior.
    """

    if not 0 <= propensity <= 1:
        raise ValueError("propensity must be in [0, 1]")
    if strength <= 0:
        raise ValueError("strength must be positive")
    m = min(max(propensity, _EPS), 1 - _EPS)
    return BetaPrior(alpha=m * strength, beta=(1 - m) * strength)


def beta_prior_from_mean_var(mean: float, variance: float) -> BetaPrior:
    """Method-of-moments Beta prior from a predicted mean and variance.

    Requires ``0 < variance < mean * (1 - mean)``. The predicted variance is
    clamped to keep the concentration positive, so an over-confident (tiny) or
    under-confident (near-maximal) variance still yields a valid prior.
    """

    if not 0 <= mean <= 1:
        raise ValueError("mean must be in [0, 1]")
    m = min(max(mean, _EPS), 1 - _EPS)
    max_var = m * (1 - m)
    v = min(max(variance, _EPS * max_var), max_var * (1 - _EPS))
    kappa = m * (1 - m) / v - 1.0
    kappa = max(kappa, _EPS)
    return BetaPrior(alpha=m * kappa, beta=(1 - m) * kappa)
