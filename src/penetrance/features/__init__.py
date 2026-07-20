"""Component 2 - mechanism-aware feature matrix and mechanism inference."""

from penetrance.features.mechanism import (
    infer_dn_gof_score,
    infer_mechanism_class,
    MECHANISM_CLASSES,
)
from penetrance.features.matrix import (
    FEATURE_COLUMNS,
    build_feature_matrix,
    build_gene_features,
)

__all__ = [
    "infer_dn_gof_score",
    "infer_mechanism_class",
    "MECHANISM_CLASSES",
    "FEATURE_COLUMNS",
    "build_feature_matrix",
    "build_gene_features",
]
