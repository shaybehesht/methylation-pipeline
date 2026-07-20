"""Mechanism inference: LoF/haploinsufficiency vs dominant-negative/gain-of-function.

This is the key novel feature. Rather than trusting a hand-annotated mechanism
label (which would make the downstream model circular), we *infer* a continuous
``dn_gof_score in [0, 1]`` from mechanism-revealing signals that are themselves
observable from ClinVar / gnomAD:

* ``mis_lof_path_ratio`` - fraction of pathogenic variants that are missense
  rather than truncating. Dominant-negative / gain-of-function genes are
  dominated by clustered missense variants (tubulins ~1.0); haploinsufficient
  genes carry many truncating pathogenic alleles (BRCA ~0.3).
* ``mis_oe`` - regional/missense observed-over-expected. Strong missense
  constraint (low o/e) points to a specific missense-driven mechanism.
* ``phaplo`` / ``pli`` - high haploinsufficiency probability pulls toward the
  LoF/HI class.

The weights are fixed, interpretable and deliberately dominated by the
missense-vs-LoF pathogenic ratio, which is the cleanest mechanistic signal. The
score is used as a model feature; the design's "open decision" (hard input vs
latent) is resolved by feeding this soft, learned-from-data score rather than a
one-hot oracle label.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
import pandas as pd

MECHANISM_CLASSES = ("LoF_HI", "mixed", "DN_GoF")

# Fixed logistic weights (see module docstring). Dominated by the pathogenic
# missense:LoF ratio, which is the least ambiguous mechanistic signal.
_W_RATIO = 5.0
_RATIO_CENTER = 0.55
_W_MIS_OE = 1.5
_MIS_OE_CENTER = 0.5
_W_PHAPLO = 1.5
_PHAPLO_CENTER = 0.55


def _sigmoid(z: np.ndarray | float):
    return 1.0 / (1.0 + np.exp(-z))


def infer_dn_gof_score(features: Mapping[str, float] | pd.DataFrame):
    """Continuous score in [0, 1]; 1 = dominant-negative/GoF, 0 = LoF/HI.

    Accepts a single mapping (returns a float) or a DataFrame (returns a Series).
    """

    if isinstance(features, pd.DataFrame):
        ratio = features["mis_lof_path_ratio"].to_numpy(dtype=float)
        mis_oe = features["mis_oe"].to_numpy(dtype=float)
        phaplo = features["phaplo"].to_numpy(dtype=float)
        z = (
            _W_RATIO * (ratio - _RATIO_CENTER)
            + _W_MIS_OE * (_MIS_OE_CENTER - mis_oe)
            - _W_PHAPLO * (phaplo - _PHAPLO_CENTER)
        )
        return pd.Series(_sigmoid(z), index=features.index, name="dn_gof_score")

    z = (
        _W_RATIO * (float(features["mis_lof_path_ratio"]) - _RATIO_CENTER)
        + _W_MIS_OE * (_MIS_OE_CENTER - float(features["mis_oe"]))
        - _W_PHAPLO * (float(features["phaplo"]) - _PHAPLO_CENTER)
    )
    return float(_sigmoid(z))


def infer_mechanism_class(features: Mapping[str, float] | pd.DataFrame):
    """Discretise the DN/GoF score into ``MECHANISM_CLASSES``."""

    score = infer_dn_gof_score(features)

    def _classify(s: float) -> str:
        if s >= 0.66:
            return "DN_GoF"
        if s <= 0.34:
            return "LoF_HI"
        return "mixed"

    if isinstance(score, pd.Series):
        return score.map(_classify)
    return _classify(score)
