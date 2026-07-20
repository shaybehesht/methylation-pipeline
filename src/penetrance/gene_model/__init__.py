"""Component 3 - the learned, mechanism-aware gene-level penetrance model.

This is the package's core contribution: a calibrated regression that predicts a
gene's *incomplete-penetrance propensity* (a continuous penetrance value in
``[0, 1]``) from mechanism, constraint and the pathogenic allele-frequency
spectrum, and then maps that prediction to a Beta prior ``(alpha, beta)`` for the
per-variant layer. Gene-family-aware cross-validation prevents paralog leakage
(e.g. tubulins co-occurring in train and test).
"""

from penetrance.gene_model.prior import (
    BetaPrior,
    propensity_to_beta_prior,
    beta_prior_from_mean_var,
)
from penetrance.gene_model.cv import GeneFamilyKFold, family_aware_split
from penetrance.gene_model.model import (
    GenePenetranceModel,
    train_gene_model,
)

__all__ = [
    "BetaPrior",
    "propensity_to_beta_prior",
    "beta_prior_from_mean_var",
    "GeneFamilyKFold",
    "family_aware_split",
    "GenePenetranceModel",
    "train_gene_model",
]
