"""Component 4 - the off-the-shelf per-variant Bayesian layer.

None of the estimators here are novel: the Beta-Binomial posterior is CalPen /
Kroncke (PLoS Genet 2020), the case/control + prevalence formula is the
ADpenetrance model (Wright/KCL), and the maximum-credible-allele-frequency bound
is Whiffin/Ware (Genet Med 2017). The package reuses them and, crucially, lets
the learned gene model (Component 3) supply the Beta prior ``(alpha, beta)``
instead of a flat one.
"""

from penetrance.variant.betabinom import (
    BetaBinomialPosterior,
    beta_binomial_posterior,
)
from penetrance.variant.casecontrol import (
    case_control_penetrance,
    case_control_penetrance_from_counts,
)
from penetrance.variant.maxaf import (
    max_credible_af,
    penetrance_upper_bound_from_af,
)
from penetrance.variant.estimator import (
    VariantPenetranceEstimate,
    estimate_variant_penetrance,
)

__all__ = [
    "BetaBinomialPosterior",
    "beta_binomial_posterior",
    "case_control_penetrance",
    "case_control_penetrance_from_counts",
    "max_credible_af",
    "penetrance_upper_bound_from_af",
    "VariantPenetranceEstimate",
    "estimate_variant_penetrance",
]
