"""Component 5 - pluggable carrier-count sources for the per-variant layer.

Two adapters ship in v1, both implementing the same :class:`CountAdapter`
interface so they can be combined:

* :class:`FrequencyCountAdapter` - ClinVar/gnomAD-derived assembled carrier
  counts (the default, offline data source).
* :class:`LiteratureCountAdapter` - a fork-of-GeneVariantFetcher literature
  miner: it parses documents, extracts carrier counts with a *pluggable*
  extractor (LLM or the built-in rule-based fallback), gates on variant
  matching, normalises phenotype ontology terms, and keeps provenance quotes.
"""

from penetrance.adapters.base import (
    CarrierCounts,
    CountAdapter,
    combine_counts,
)
from penetrance.adapters.clinvar_gnomad import FrequencyCountAdapter
from penetrance.adapters.literature import (
    LiteratureCountAdapter,
    LiteratureRecord,
    regex_carrier_extractor,
)

__all__ = [
    "CarrierCounts",
    "CountAdapter",
    "combine_counts",
    "FrequencyCountAdapter",
    "LiteratureCountAdapter",
    "LiteratureRecord",
    "regex_carrier_extractor",
]
