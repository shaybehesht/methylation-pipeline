"""Adapter interface and carrier-count aggregation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Iterable, List, Optional


@dataclass
class CarrierCounts:
    """Affected/unaffected carrier counts for one variant, with provenance."""

    variant_id: str
    affected: float
    unaffected: float
    source: str
    allele_frequency: Optional[float] = None
    provenance: List[str] = field(default_factory=list)

    @property
    def n_carriers(self) -> float:
        return self.affected + self.unaffected

    @property
    def observed_penetrance(self) -> Optional[float]:
        n = self.n_carriers
        return self.affected / n if n > 0 else None


class CountAdapter(ABC):
    """Base class for a pluggable carrier-count source."""

    name: str = "adapter"

    @abstractmethod
    def fetch(self, variant_id: str) -> Optional[CarrierCounts]:
        """Return counts for a single variant, or ``None`` if unavailable."""

    def fetch_many(self, variant_ids: Iterable[str]) -> List[CarrierCounts]:
        out = []
        for vid in variant_ids:
            counts = self.fetch(vid)
            if counts is not None:
                out.append(counts)
        return out


def combine_counts(counts: Iterable[CarrierCounts]) -> Optional[CarrierCounts]:
    """Sum carrier counts for the same variant across adapters.

    Allele frequency is taken from the first adapter that reports one (gnomAD is
    the authoritative population frequency); provenance strings are concatenated.
    """

    counts = list(counts)
    if not counts:
        return None
    vid = counts[0].variant_id
    affected = sum(c.affected for c in counts)
    unaffected = sum(c.unaffected for c in counts)
    af = next((c.allele_frequency for c in counts if c.allele_frequency is not None), None)
    provenance: List[str] = []
    sources: List[str] = []
    for c in counts:
        sources.append(c.source)
        provenance.extend(c.provenance)
    return CarrierCounts(
        variant_id=vid,
        affected=affected,
        unaffected=unaffected,
        source="+".join(dict.fromkeys(sources)),
        allele_frequency=af,
        provenance=provenance,
    )
