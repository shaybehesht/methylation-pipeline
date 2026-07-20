"""End-to-end orchestration: gene prior -> per-variant posterior."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import pandas as pd

from penetrance.adapters.base import CarrierCounts, CountAdapter, combine_counts
from penetrance.adapters.clinvar_gnomad import FrequencyCountAdapter
from penetrance.features.matrix import build_feature_matrix, build_gene_features
from penetrance.gene_model.model import GenePenetranceModel, train_gene_model
from penetrance.gene_model.prior import BetaPrior, propensity_to_beta_prior
from penetrance.labels.loader import LabelSet, load_labels
from penetrance.variant.estimator import VariantPenetranceEstimate, estimate_variant_penetrance

# Fallback prevalence when a variant's phenotype prevalence is unknown.
_DEFAULT_PREVALENCE = 1e-3


@dataclass
class GenePrediction:
    gene: str
    propensity: float
    std: float
    prior: BetaPrior


class PenetrancePipeline:
    """Trains the gene model and produces mechanism-aware variant estimates."""

    def __init__(
        self,
        labels: Optional[LabelSet] = None,
        adapters: Optional[Sequence[CountAdapter]] = None,
        prior_strength: float = 12.0,
    ):
        self.labels = labels or load_labels()
        self.adapters: List[CountAdapter] = (
            list(adapters) if adapters is not None else [FrequencyCountAdapter(self.labels.variants)]
        )
        self.prior_strength = prior_strength
        self.model: Optional[GenePenetranceModel] = None
        self.cv_result_ = None

    # ---------------------------------------------------------------- gene
    def fit(self, **kwargs) -> "PenetrancePipeline":
        self.model, self.cv_result_ = train_gene_model(
            self.labels.genes, prior_strength=self.prior_strength, **kwargs
        )
        return self

    def _gene_row(self, gene: str) -> Optional[pd.Series]:
        return self.labels.gene_row(gene)

    def predict_gene(self, gene: str) -> Optional[GenePrediction]:
        if self.model is None:
            raise RuntimeError("call fit() before predicting")
        row = self._gene_row(gene)
        if row is None:
            return None
        X, _ = build_feature_matrix(pd.DataFrame([row]))
        mean, std = self.model.predict_with_uncertainty(X)
        prior = self.model.beta_prior_for(X)[0]
        return GenePrediction(gene=gene, propensity=float(mean[0]), std=float(std[0]), prior=prior)

    def gene_prior(self, gene: str, flat: bool = False) -> BetaPrior:
        """Beta prior for a gene: mechanism-aware (default) or flat ``Beta(1,1)``."""

        if flat:
            return BetaPrior(1.0, 1.0)
        pred = self.predict_gene(gene)
        if pred is None:
            return BetaPrior(1.0, 1.0)
        return pred.prior

    # ------------------------------------------------------------- variant
    def carrier_counts(self, variant_id: str) -> Optional[CarrierCounts]:
        collected = [c for a in self.adapters if (c := a.fetch(variant_id)) is not None]
        return combine_counts(collected)

    def estimate_variant(
        self,
        variant_id: str,
        use_gene_prior: bool = True,
        apply_af_bound: bool = True,
        override_counts: Optional[CarrierCounts] = None,
    ) -> Optional[VariantPenetranceEstimate]:
        vrow = self.labels.variants[self.labels.variants["variant_id"] == variant_id]
        if vrow.empty:
            return None
        vrow = vrow.iloc[0]
        gene = vrow["gene"]

        counts = override_counts or self.carrier_counts(variant_id)
        if counts is None:
            counts = CarrierCounts(variant_id, 0.0, 0.0, "none")

        prior = self.gene_prior(gene, flat=not use_gene_prior)

        af = counts.allele_frequency
        if af is None and "gnomad_af" in vrow:
            af = float(vrow["gnomad_af"]) if not pd.isna(vrow["gnomad_af"]) else None
        prevalence = float(vrow["prevalence_K"]) if "prevalence_K" in vrow and not pd.isna(vrow["prevalence_K"]) else _DEFAULT_PREVALENCE

        return estimate_variant_penetrance(
            affected=counts.affected,
            unaffected=counts.unaffected,
            prior_alpha=prior.alpha,
            prior_beta=prior.beta,
            allele_frequency=af,
            prevalence=prevalence,
            apply_af_bound=apply_af_bound,
        )
