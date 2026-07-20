"""Mechanism-aware penetrance prediction.

The package is organised around the architecture described in the design doc:

* :mod:`penetrance.labels` - curated, literature-derived penetrance ground truth
  (Component 1).
* :mod:`penetrance.features` - the mechanism-aware feature matrix and mechanism
  inference (Component 2).
* :mod:`penetrance.gene_model` - the learned gene-level incomplete-penetrance
  propensity model and its mapping to a Beta prior ``(alpha, beta)``
  (Component 3). This is the core contribution.
* :mod:`penetrance.variant` - the off-the-shelf per-variant Bayesian layer:
  Beta-Binomial posterior, case/control-prevalence Bayes and the Whiffin/Ware
  maximum-credible-allele-frequency bound (Component 4).
* :mod:`penetrance.adapters` - pluggable carrier-count sources (Component 5).

The high-level orchestration lives in :mod:`penetrance.pipeline`.
"""

from __future__ import annotations

__version__ = "0.1.0"

from penetrance.variant.betabinom import (
    BetaBinomialPosterior,
    beta_binomial_posterior,
)
from penetrance.variant.casecontrol import case_control_penetrance
from penetrance.variant.maxaf import max_credible_af, penetrance_upper_bound_from_af
from penetrance.gene_model.prior import propensity_to_beta_prior, BetaPrior

__all__ = [
    "__version__",
    "BetaBinomialPosterior",
    "beta_binomial_posterior",
    "case_control_penetrance",
    "max_credible_af",
    "penetrance_upper_bound_from_af",
    "propensity_to_beta_prior",
    "BetaPrior",
]
