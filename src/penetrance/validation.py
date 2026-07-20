"""Component 6 - validation that proves the point.

Three questions, matching the design:

1. Does the gene model separate high-penetrance dominant-negative/GoF genes
   (tubulins) from incomplete-penetrance haploinsufficiency genes
   (BRCA/LDLR/SDHx)? -> :func:`mechanism_separation`.
2. Does the mechanism-aware Beta prior beat a flat prior for *sparse* variants
   (few carriers)? -> :func:`sparse_prior_benefit`, a mask-and-recover
   experiment.
3. Are the gene model's cross-validated predictions calibrated? ->
   :func:`gene_model_calibration`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from penetrance.gene_model.prior import BetaPrior
from penetrance.pipeline import PenetrancePipeline
from penetrance.variant.betabinom import beta_binomial_posterior


def mechanism_separation(
    pipeline: PenetrancePipeline,
    high_penetrance_families: Sequence[str] = ("tubulin", "collagen", "fibrillin"),
    incomplete_families: Sequence[str] = ("hboc", "fh", "pgl_pheo"),
) -> Dict[str, float]:
    """Compare predicted propensity for high- vs incomplete-penetrance families."""

    genes = pipeline.labels.genes
    preds = {g: pipeline.predict_gene(g) for g in genes["gene"]}

    def _mean_for(families):
        vals = [
            preds[g].propensity
            for g, fam in zip(genes["gene"], genes["gene_family"])
            if fam in families and preds[g] is not None
        ]
        return float(np.mean(vals)) if vals else float("nan")

    high = _mean_for(set(high_penetrance_families))
    low = _mean_for(set(incomplete_families))
    return {
        "high_penetrance_mean_propensity": high,
        "incomplete_penetrance_mean_propensity": low,
        "separation": high - low,
    }


@dataclass
class SparsePriorResult:
    table: pd.DataFrame
    summary: pd.DataFrame

    def improvement(self) -> pd.Series:
        """Fractional MAE reduction of the mechanism prior over the flat prior."""

        s = self.summary
        return (s["flat_mae"] - s["mechanism_mae"]) / s["flat_mae"]


def sparse_prior_benefit(
    pipeline: PenetrancePipeline,
    subsample_sizes: Sequence[int] = (1, 2, 3, 5, 8),
    n_trials: int = 400,
    min_carriers: int = 30,
    random_state: int = 0,
) -> SparsePriorResult:
    """Mask-and-recover experiment on sparse variants.

    For each well-observed variant we treat its empirical penetrance as truth,
    draw ``k`` carriers without replacement to mimic a sparse variant, and
    estimate penetrance with (a) a flat ``Beta(1,1)`` prior and (b) the
    mechanism-aware gene prior. Lower error for (b) at small ``k`` is the whole
    justification for the gene model.
    """

    rng = np.random.default_rng(random_state)
    variants = pipeline.labels.variants
    rows: List[dict] = []

    for _, v in variants.iterrows():
        n = int(v["affected_carriers"] + v["unaffected_carriers"])
        if n < min_carriers:
            continue
        a_full = int(v["affected_carriers"])
        p_true = a_full / n
        pool = np.array([1] * a_full + [0] * (n - a_full))
        prior = pipeline.gene_prior(v["gene"], flat=False)
        flat = BetaPrior(1.0, 1.0)

        for k in subsample_sizes:
            if k > n:
                continue
            for _ in range(n_trials):
                draw = rng.choice(pool, size=k, replace=False)
                a_k = int(draw.sum())
                u_k = k - a_k
                est_flat = beta_binomial_posterior(a_k, u_k, flat.alpha, flat.beta).mean
                est_mech = beta_binomial_posterior(a_k, u_k, prior.alpha, prior.beta).mean
                rows.append(
                    {
                        "variant_id": v["variant_id"],
                        "gene": v["gene"],
                        "k": k,
                        "p_true": p_true,
                        "err_flat": abs(est_flat - p_true),
                        "err_mech": abs(est_mech - p_true),
                    }
                )

    table = pd.DataFrame(rows)
    summary = (
        table.groupby("k")
        .agg(flat_mae=("err_flat", "mean"), mechanism_mae=("err_mech", "mean"))
        .reset_index()
    )
    return SparsePriorResult(table=table, summary=summary)


def gene_model_calibration(
    pipeline: PenetrancePipeline, n_bins: int = 5
) -> Dict[str, object]:
    """Calibration of cross-validated gene predictions vs true penetrance."""

    cv = pipeline.cv_result_
    if cv is None:
        raise RuntimeError("pipeline has no CV result; call fit() first")
    y_true = np.asarray(cv.y_true)
    y_pred = np.asarray(cv.y_pred)
    order = np.argsort(y_pred)
    bins = np.array_split(order, n_bins)
    curve = []
    for b in bins:
        if len(b) == 0:
            continue
        curve.append(
            {
                "pred_mean": float(y_pred[b].mean()),
                "true_mean": float(y_true[b].mean()),
                "n": int(len(b)),
            }
        )
    curve_df = pd.DataFrame(curve)
    ece = float(np.average((curve_df["pred_mean"] - curve_df["true_mean"]).abs(), weights=curve_df["n"]))
    return {"metrics": cv.metrics, "reliability": curve_df, "expected_calibration_error": ece}
