"""Component 1 - curated penetrance ground truth.

The gene-level and variant-level penetrance labels are stored as packaged CSVs
under ``penetrance/data``. They are curated approximations drawn from the
population/biobank and clinical literature cited in each row's ``source`` field
(Wright et al. Nat Genet 2024; Forrest/Huang et al. Nat Genet 2025; the Science
2024/2025 ML-penetrance gene set; ClinGen/OMIM; and the individual studies named
per gene). They exist so the pipeline is runnable offline and reproducible; in
production the same tables are meant to be regenerated from primary sources.
"""

from penetrance.labels.loader import (
    LabelSet,
    load_gene_labels,
    load_variant_labels,
    load_labels,
    label_weight,
)

__all__ = [
    "LabelSet",
    "load_gene_labels",
    "load_variant_labels",
    "load_labels",
    "label_weight",
]
