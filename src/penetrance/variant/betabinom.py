"""Beta-Binomial posterior for per-variant penetrance.

Penetrance is ``p = P(affected | carrier)``. Given ``a`` affected and ``u``
unaffected carriers and a ``Beta(alpha, beta)`` prior, the posterior is

    p | data ~ Beta(alpha + a, beta + u)

with posterior mean ``(alpha + a) / (alpha + beta + a + u)``. This is the exact
conjugate update used by CalPen / Kroncke (PLoS Genet 2020). The only thing this
package changes is *where the prior comes from*: the mechanism-aware gene model
supplies ``(alpha, beta)`` (see :mod:`penetrance.gene_model.prior`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

from scipy import stats


@dataclass(frozen=True)
class BetaBinomialPosterior:
    """A Beta posterior over a variant's penetrance.

    Attributes
    ----------
    alpha, beta:
        Posterior Beta parameters (already updated with the observed counts).
    prior_alpha, prior_beta:
        The prior parameters, retained for provenance / diagnostics.
    affected, unaffected:
        The observed carrier counts that produced the posterior.
    """

    alpha: float
    beta: float
    prior_alpha: float
    prior_beta: float
    affected: float
    unaffected: float

    @property
    def mean(self) -> float:
        return self.alpha / (self.alpha + self.beta)

    @property
    def mode(self) -> float:
        # Mode is only defined for alpha, beta > 1; fall back to the mean.
        if self.alpha > 1 and self.beta > 1:
            return (self.alpha - 1) / (self.alpha + self.beta - 2)
        return self.mean

    @property
    def variance(self) -> float:
        a, b = self.alpha, self.beta
        return (a * b) / ((a + b) ** 2 * (a + b + 1))

    @property
    def n_observed(self) -> float:
        return self.affected + self.unaffected

    @property
    def prior_strength(self) -> float:
        """Effective sample size of the prior (pseudo-count)."""

        return self.prior_alpha + self.prior_beta

    def credible_interval(self, level: float = 0.95) -> Tuple[float, float]:
        """Equal-tailed Bayesian credible interval at the requested level."""

        if not 0 < level < 1:
            raise ValueError("level must be in (0, 1)")
        tail = (1.0 - level) / 2.0
        dist = stats.beta(self.alpha, self.beta)
        return float(dist.ppf(tail)), float(dist.ppf(1.0 - tail))

    def prob_penetrance_above(self, threshold: float) -> float:
        """Posterior probability that penetrance exceeds ``threshold``."""

        return float(stats.beta(self.alpha, self.beta).sf(threshold))


def beta_binomial_posterior(
    affected: float,
    unaffected: float,
    prior_alpha: float = 1.0,
    prior_beta: float = 1.0,
) -> BetaBinomialPosterior:
    """Conjugate Beta-Binomial update.

    Parameters
    ----------
    affected, unaffected:
        Carrier counts. Ideally genotype-first ascertained to avoid inflating
        the affected count.
    prior_alpha, prior_beta:
        Beta prior parameters. Defaults to the flat ``Beta(1, 1)`` prior; pass
        the mechanism-aware gene prior for sparse variants.
    """

    if affected < 0 or unaffected < 0:
        raise ValueError("carrier counts must be non-negative")
    if prior_alpha <= 0 or prior_beta <= 0:
        raise ValueError("prior parameters must be positive")
    return BetaBinomialPosterior(
        alpha=prior_alpha + affected,
        beta=prior_beta + unaffected,
        prior_alpha=prior_alpha,
        prior_beta=prior_beta,
        affected=affected,
        unaffected=unaffected,
    )
