"""Assemble the mechanism-aware feature matrix used by the gene model."""

from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd

from penetrance.features.mechanism import infer_dn_gof_score

# The ordered feature set fed to the gene model. Grouped by the three pillars of
# the thesis: mechanism, constraint, and the pathogenic allele-frequency
# spectrum, plus modifier/expressivity proxies.
FEATURE_COLUMNS: List[str] = [
    # mechanism
    "dn_gof_score",
    "mis_lof_path_ratio",
    # constraint
    "loeuf",
    "pli",
    "shet",
    "phaplo",
    "ptriplo",
    "mis_oe",
    "haplo_triplo_asymmetry",
    # pathogenic allele-frequency spectrum
    "log_path_af_mean",
    "log_path_af_max",
    "af_spectrum_tail",
    "log_n_path",
    # modifiers / expressivity
    "hpo_breadth",
    "log_paralog_count",
]

_LOG_FLOOR = 1e-9


def build_gene_features(genes: pd.DataFrame) -> pd.DataFrame:
    """Derive the engineered feature columns from the raw gene table.

    The input must contain the raw columns present in ``data/genes.csv`` (see
    :data:`penetrance.labels.loader.GENE_FEATURE_COLUMNS`). Returns a new frame
    indexed like the input with ``gene``, ``gene_family`` and all
    :data:`FEATURE_COLUMNS`.
    """

    df = genes.copy()
    df["dn_gof_score"] = infer_dn_gof_score(df)
    df["haplo_triplo_asymmetry"] = df["phaplo"] - df["ptriplo"]
    df["log_path_af_mean"] = np.log10(df["path_af_mean"].clip(lower=_LOG_FLOOR))
    df["log_path_af_max"] = np.log10(df["path_af_max"].clip(lower=_LOG_FLOOR))
    # "fat tail" signature: how much higher the most-common pathogenic allele is
    # than the average one. A large gap = a too-common pathogenic allele exists.
    df["af_spectrum_tail"] = df["log_path_af_max"] - df["log_path_af_mean"]
    df["log_n_path"] = np.log10(df["n_path"].clip(lower=1))
    df["log_paralog_count"] = np.log10(df["paralog_count"].clip(lower=1) + 1)

    keep = ["gene", "gene_family"] + FEATURE_COLUMNS
    return df[keep].reset_index(drop=True)


def build_feature_matrix(genes: pd.DataFrame):
    """Return ``(X, meta)`` where ``X`` is the model matrix and ``meta`` has ids.

    ``meta`` carries ``gene`` and ``gene_family`` (the latter is the grouping key
    for family-aware CV).
    """

    feats = build_gene_features(genes)
    X = feats[FEATURE_COLUMNS].astype(float)
    meta = feats[["gene", "gene_family"]].copy()
    return X, meta
