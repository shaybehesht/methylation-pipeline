"""Loaders for the curated penetrance label tables."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from typing import Optional

import pandas as pd

# Ascertainment / source-type -> label variance. Population/biobank estimates
# are the least ascertainment-biased and therefore the most trustworthy targets;
# single-report clinical series are the noisiest. Variance is expressed on the
# penetrance scale (units of probability^2) and is turned into a regression
# sample weight of 1 / variance.
SOURCE_TYPE_VARIANCE = {
    "population": 0.01,
    "clinical": 0.03,
    "family": 0.05,
    "single_study": 0.08,
}
_DEFAULT_VARIANCE = 0.05

GENE_FEATURE_COLUMNS = [
    "loeuf",
    "pli",
    "shet",
    "phaplo",
    "ptriplo",
    "mis_oe",
    "mis_lof_path_ratio",
    "path_af_mean",
    "path_af_max",
    "n_path",
    "hpo_breadth",
    "paralog_count",
]


@dataclass
class LabelSet:
    """Container bundling the gene-level and variant-level label tables."""

    genes: pd.DataFrame
    variants: pd.DataFrame

    def gene_row(self, gene: str) -> Optional[pd.Series]:
        hit = self.genes[self.genes["gene"] == gene]
        if hit.empty:
            return None
        return hit.iloc[0]


def _data_path(name: str):
    return resources.files("penetrance.data").joinpath(name)


def label_weight(source_type: str) -> float:
    """Regression sample weight for a label given its source/ascertainment type.

    Noisier labels (higher variance) contribute less. Implements the locked
    decision to "use the midpoint with a per-label weight/variance so noisy
    labels count less".
    """

    variance = SOURCE_TYPE_VARIANCE.get(source_type, _DEFAULT_VARIANCE)
    return 1.0 / variance


def load_gene_labels() -> pd.DataFrame:
    """Load the gene-level penetrance labels + curated mechanism/constraint features."""

    with resources.as_file(_data_path("genes.csv")) as path:
        df = pd.read_csv(path)
    df["penetrance_variance"] = df["penetrance_source_type"].map(
        lambda s: SOURCE_TYPE_VARIANCE.get(s, _DEFAULT_VARIANCE)
    )
    df["label_weight"] = df["penetrance_source_type"].map(label_weight)
    return df


def load_variant_labels() -> pd.DataFrame:
    """Load variant-level carrier counts and observed penetrance."""

    with resources.as_file(_data_path("variants.csv")) as path:
        df = pd.read_csv(path)
    total = df["affected_carriers"] + df["unaffected_carriers"]
    df["n_carriers"] = total
    df["penetrance_obs"] = df["affected_carriers"] / total.where(total > 0, other=pd.NA)
    return df


def load_labels() -> LabelSet:
    """Load both label tables into a :class:`LabelSet`."""

    return LabelSet(genes=load_gene_labels(), variants=load_variant_labels())
