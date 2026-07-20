"""ClinVar/gnomAD frequency-based carrier-count adapter.

This adapter treats an assembled frequency/count table as its backing store. By
default it loads the packaged ``data/variants.csv`` (curated ClinVar/gnomAD-style
counts) so the pipeline runs offline; point it at a table pulled from live
ClinVar + gnomAD to use real data. The table must expose ``variant_id``,
``affected_carriers``, ``unaffected_carriers`` and (optionally) ``gnomad_af``.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from penetrance.adapters.base import CarrierCounts, CountAdapter
from penetrance.labels.loader import load_variant_labels


class FrequencyCountAdapter(CountAdapter):
    name = "clinvar_gnomad"

    def __init__(self, table: Optional[pd.DataFrame] = None):
        if table is None:
            table = load_variant_labels()
        required = {"variant_id", "affected_carriers", "unaffected_carriers"}
        missing = required - set(table.columns)
        if missing:
            raise ValueError(f"frequency table is missing columns: {sorted(missing)}")
        self._table = table.set_index("variant_id", drop=False)

    def fetch(self, variant_id: str) -> Optional[CarrierCounts]:
        if variant_id not in self._table.index:
            return None
        row = self._table.loc[variant_id]
        if isinstance(row, pd.DataFrame):  # duplicate ids -> take the first
            row = row.iloc[0]
        af = row.get("gnomad_af") if hasattr(row, "get") else None
        return CarrierCounts(
            variant_id=variant_id,
            affected=float(row["affected_carriers"]),
            unaffected=float(row["unaffected_carriers"]),
            source=self.name,
            allele_frequency=(float(af) if af is not None and not pd.isna(af) else None),
            provenance=[f"{self.name}: {row.get('source', 'n/a')}"],
        )
