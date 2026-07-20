"""The learned gene-level incomplete-penetrance-propensity model.

A calibrated gradient-boosted regression with a *continuous* penetrance target.
It exposes:

* out-of-fold predictions from gene-family-aware CV (for honest validation and
  calibration),
* per-gene predictive uncertainty from quantile regressors,
* SHAP attributions (optional), and
* a direct mapping to a Beta prior for the per-variant layer.

LightGBM is used when available; otherwise it falls back to scikit-learn's
``HistGradientBoostingRegressor`` so the package works with no extra deps.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

from penetrance.features.matrix import FEATURE_COLUMNS, build_feature_matrix
from penetrance.gene_model.cv import GeneFamilyKFold
from penetrance.gene_model.prior import (
    BetaPrior,
    beta_prior_from_mean_var,
    propensity_to_beta_prior,
)

try:  # pragma: no cover - exercised indirectly
    import lightgbm as lgb

    _HAS_LGB = True
except Exception:  # pragma: no cover
    _HAS_LGB = False

from sklearn.ensemble import HistGradientBoostingRegressor


# Small-data-friendly hyper-parameters. The label set is ~tens of genes, so we
# keep trees shallow and heavily regularised to learn mechanism, not memorise.
_LGB_PARAMS = dict(
    objective="regression",
    n_estimators=300,
    learning_rate=0.03,
    num_leaves=7,
    min_child_samples=4,
    subsample=0.8,
    subsample_freq=1,
    colsample_bytree=0.8,
    reg_lambda=1.0,
    verbose=-1,
)
_SK_PARAMS = dict(
    max_iter=300,
    learning_rate=0.03,
    max_leaf_nodes=7,
    min_samples_leaf=4,
    l2_regularization=1.0,
)


def _make_regressor(quantile: Optional[float] = None):
    if _HAS_LGB:
        params = dict(_LGB_PARAMS)
        if quantile is not None:
            params["objective"] = "quantile"
            params["alpha"] = quantile
        return lgb.LGBMRegressor(**params)
    params = dict(_SK_PARAMS)
    if quantile is not None:
        params["loss"] = "quantile"
        params["quantile"] = quantile
    else:
        params["loss"] = "squared_error"
    return HistGradientBoostingRegressor(**params)


@dataclass
class GeneModelCVResult:
    genes: List[str]
    families: List[str]
    y_true: np.ndarray
    y_pred: np.ndarray
    weights: np.ndarray
    metrics: Dict[str, float] = field(default_factory=dict)


class GenePenetranceModel:
    """Mechanism-aware gene-level penetrance-propensity regressor."""

    def __init__(
        self,
        prior_strength: float = 12.0,
        estimate_uncertainty: bool = True,
        calibrate: bool = True,
        n_splits: int = 5,
        random_state: int = 0,
    ):
        self.prior_strength = prior_strength
        self.estimate_uncertainty = estimate_uncertainty
        self.calibrate = calibrate
        self.n_splits = n_splits
        self.random_state = random_state

        self.estimator_ = None
        self.lower_estimator_ = None
        self.upper_estimator_ = None
        self.calibrator_: Optional[IsotonicRegression] = None
        self.feature_columns_: List[str] = list(FEATURE_COLUMNS)
        self.cv_result_: Optional[GeneModelCVResult] = None
        self.backend_ = "lightgbm" if _HAS_LGB else "sklearn"

    # ------------------------------------------------------------------ fit
    def fit(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        sample_weight: Optional[np.ndarray] = None,
        groups: Optional[np.ndarray] = None,
        run_cv: bool = True,
    ) -> "GenePenetranceModel":
        X = X[self.feature_columns_].astype(float)
        y = np.asarray(y, dtype=float)
        if sample_weight is None:
            sample_weight = np.ones(len(y))
        sample_weight = np.asarray(sample_weight, dtype=float)

        if run_cv and groups is not None:
            self.cv_result_ = self._run_cv(X, y, sample_weight, groups)
            if self.calibrate:
                self._fit_calibrator(
                    self.cv_result_.y_pred, self.cv_result_.y_true, self.cv_result_.weights
                )

        self.estimator_ = _make_regressor()
        self.estimator_.fit(X, y, sample_weight=sample_weight)
        if self.estimate_uncertainty:
            self.lower_estimator_ = _make_regressor(quantile=0.16)
            self.upper_estimator_ = _make_regressor(quantile=0.84)
            self.lower_estimator_.fit(X, y, sample_weight=sample_weight)
            self.upper_estimator_.fit(X, y, sample_weight=sample_weight)
        return self

    def _fit_calibrator(self, oof_pred, y_true, weights) -> None:
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        iso.fit(oof_pred, y_true, sample_weight=weights)
        self.calibrator_ = iso

    def _run_cv(self, X, y, sample_weight, groups) -> GeneModelCVResult:
        splitter = GeneFamilyKFold(
            n_splits=self.n_splits, shuffle=True, random_state=self.random_state
        )
        oof = np.full(len(y), np.nan)
        for train_idx, test_idx in splitter.split(X, y, groups):
            reg = _make_regressor()
            reg.fit(
                X.iloc[train_idx], y[train_idx], sample_weight=sample_weight[train_idx]
            )
            oof[test_idx] = reg.predict(X.iloc[test_idx])
        mask = ~np.isnan(oof)
        oof = np.clip(oof, 0.0, 1.0)
        result = GeneModelCVResult(
            genes=list(np.asarray(groups)[mask]),
            families=list(np.asarray(groups)[mask]),
            y_true=y[mask],
            y_pred=oof[mask],
            weights=sample_weight[mask],
            metrics=_regression_metrics(y[mask], oof[mask], sample_weight[mask]),
        )
        return result

    # -------------------------------------------------------------- predict
    def _apply_calibration(self, raw: np.ndarray) -> np.ndarray:
        if self.calibrator_ is not None:
            return np.clip(self.calibrator_.predict(raw), 0.0, 1.0)
        return np.clip(raw, 0.0, 1.0)

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        if self.estimator_ is None:
            raise RuntimeError("model is not fitted")
        X = X[self.feature_columns_].astype(float)
        return self._apply_calibration(self.estimator_.predict(X))

    def predict_with_uncertainty(self, X: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Return calibrated ``(mean, std)`` per gene.

        Std comes from the 16/84 quantile spread (a Gaussian-equivalent sigma).
        When quantile models are disabled, a small floor std is returned so the
        Beta prior stays proper.
        """

        mean = self.predict(X)
        if not self.estimate_uncertainty or self.lower_estimator_ is None:
            return mean, np.full_like(mean, 0.15)
        Xf = X[self.feature_columns_].astype(float)
        lo = np.clip(self.lower_estimator_.predict(Xf), 0.0, 1.0)
        hi = np.clip(self.upper_estimator_.predict(Xf), 0.0, 1.0)
        std = np.clip((hi - lo) / 2.0, 0.02, 0.45)
        return mean, std

    # ------------------------------------------------------------ to prior
    def beta_prior_for(
        self, X: pd.DataFrame, use_uncertainty: bool = False
    ) -> List[BetaPrior]:
        """Beta prior(s) for the per-variant layer, one per row of ``X``.

        By default the prior uses a fixed concentration (``prior_strength``),
        which is stable on small label sets. Set ``use_uncertainty=True`` to
        instead derive the concentration from the model's per-gene predictive
        variance (method-of-moments); this is noisier with few training genes.
        """

        if use_uncertainty and self.estimate_uncertainty:
            mean, std = self.predict_with_uncertainty(X)
            return [beta_prior_from_mean_var(m, s ** 2) for m, s in zip(mean, std)]
        mean = self.predict(X)
        return [propensity_to_beta_prior(m, self.prior_strength) for m in mean]

    # ------------------------------------------------------------- explain
    def shap_values(self, X: pd.DataFrame):
        """SHAP values for the (uncalibrated) mean model. Requires ``shap``."""

        import shap  # local import; optional dependency

        Xf = X[self.feature_columns_].astype(float)
        explainer = shap.TreeExplainer(self.estimator_)
        return explainer.shap_values(Xf)

    def feature_importance(self) -> pd.Series:
        if self.estimator_ is None:
            raise RuntimeError("model is not fitted")
        if hasattr(self.estimator_, "feature_importances_"):
            imp = np.asarray(self.estimator_.feature_importances_, dtype=float)
        else:  # pragma: no cover
            imp = np.zeros(len(self.feature_columns_))
        s = pd.Series(imp, index=self.feature_columns_).sort_values(ascending=False)
        return s


def _regression_metrics(y_true, y_pred, weights) -> Dict[str, float]:
    err = y_pred - y_true
    mae = float(np.average(np.abs(err), weights=weights))
    rmse = float(np.sqrt(np.average(err ** 2, weights=weights)))
    # Weighted Pearson correlation.
    wm = lambda v: np.average(v, weights=weights)
    yt, yp = y_true - wm(y_true), y_pred - wm(y_pred)
    denom = np.sqrt(wm(yt ** 2) * wm(yp ** 2))
    corr = float(wm(yt * yp) / denom) if denom > 0 else 0.0
    # Spearman (rank) correlation, unweighted.
    from scipy.stats import spearmanr

    rho = float(spearmanr(y_true, y_pred).correlation) if len(y_true) > 2 else 0.0
    return {"mae": mae, "rmse": rmse, "pearson": corr, "spearman": rho}


def train_gene_model(
    genes: pd.DataFrame,
    prior_strength: float = 12.0,
    n_splits: int = 5,
    random_state: int = 0,
    estimate_uncertainty: bool = True,
    calibrate: bool = True,
) -> Tuple[GenePenetranceModel, GeneModelCVResult]:
    """Convenience: build features, run family-aware CV, and fit on all data.

    Returns the fitted model and its cross-validation result.
    """

    X, meta = build_feature_matrix(genes)
    y = genes["penetrance"].to_numpy(dtype=float)
    weights = genes["label_weight"].to_numpy(dtype=float) if "label_weight" in genes else None
    groups = meta["gene_family"].to_numpy()

    model = GenePenetranceModel(
        prior_strength=prior_strength,
        n_splits=n_splits,
        random_state=random_state,
        estimate_uncertainty=estimate_uncertainty,
        calibrate=calibrate,
    )
    model.fit(X, y, sample_weight=weights, groups=groups, run_cv=True)
    return model, model.cv_result_
