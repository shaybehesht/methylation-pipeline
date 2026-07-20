"""Literature-mined carrier-count adapter (fork-of-GeneVariantFetcher).

Given a corpus of documents, this adapter extracts affected/unaffected carrier
counts for a target variant. It implements the pieces the design calls for:

1. **Pluggable extraction** - an ``extractor`` callable does the parsing. Plug in
   an LLM (a function that calls your model) for production; the built-in
   :func:`regex_carrier_extractor` is a dependency-free fallback so the adapter
   runs offline and is unit-testable.
2. **Variant-matching gate** - extractions are only accepted if the document
   mentions the target variant, matched via normalised HGVS / protein tokens
   (3-letter and 1-letter amino-acid forms, cDNA change). This prevents pulling
   counts for the wrong variant.
3. **Ontology normalisation** - reported phenotype strings are mapped to a
   canonical term via a supplied dictionary.
4. **Provenance** - every accepted count keeps the sentence it came from.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import pandas as pd

from penetrance.adapters.base import CarrierCounts, CountAdapter
from penetrance.labels.loader import load_variant_labels

_AA3_TO_1 = {
    "ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C", "gln": "Q",
    "glu": "E", "gly": "G", "his": "H", "ile": "I", "leu": "L", "lys": "K",
    "met": "M", "phe": "F", "pro": "P", "ser": "S", "thr": "T", "trp": "W",
    "tyr": "Y", "val": "V", "ter": "X", "*": "X",
}

_PROT_RE = re.compile(
    r"p\.?\(?\s*([A-Za-z]{3}|[A-Z])\s*(\d+)\s*([A-Za-z]{3}|[A-Z*]|fs|Ter|del|dup)",
    re.IGNORECASE,
)
_CDNA_RE = re.compile(r"c\.\s*([0-9_+\-]+[ACGT]*>?[ACGT]*(?:del|dup|ins)?[ACGT]*)", re.IGNORECASE)


@dataclass
class LiteratureRecord:
    """A single document (abstract / full text / supplement chunk)."""

    doc_id: str
    gene: str
    text: str
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class ExtractionResult:
    variant_token: str
    affected: float
    unaffected: float
    phenotype: Optional[str]
    quote: str


# An extractor takes a record + the set of acceptable variant tokens and returns
# any counts it can find. LLM extractors have the same signature.
Extractor = Callable[[LiteratureRecord, Sequence[str]], List[ExtractionResult]]


def _norm_protein(match: re.Match) -> List[str]:
    ref, pos, alt = match.group(1), match.group(2), match.group(3)

    def to1(tok: str) -> str:
        t = tok.lower()
        if t in _AA3_TO_1:
            return _AA3_TO_1[t]
        return tok.upper()

    ref1, alt1 = to1(ref), to1(alt)
    tokens = {f"p.{ref1}{pos}{alt1}".lower()}
    tokens.add(f"{ref1}{pos}{alt1}".lower())
    return list(tokens)


def normalize_variant_tokens(hgvs: str = "", variant_id: str = "") -> List[str]:
    """Normalised match tokens for a variant (protein + cDNA forms)."""

    text = f"{hgvs} {variant_id}"
    tokens = set()
    for m in _PROT_RE.finditer(text):
        tokens.update(_norm_protein(m))
    for m in _CDNA_RE.finditer(text):
        tokens.add(f"c.{m.group(1).lower().replace(' ', '')}")
    return sorted(tokens)


def _variant_mentioned(text: str, tokens: Sequence[str]) -> Optional[str]:
    low = text.lower().replace(" ", "")
    for tok in tokens:
        if tok and tok.replace(" ", "") in low:
            return tok
    # Also try protein tokens re-extracted from the document itself.
    for m in _PROT_RE.finditer(text):
        for t in _norm_protein(m):
            if t in tokens:
                return t
    return None


_COUNT_PATTERNS = [
    # "12 of 20 carriers were affected"
    re.compile(r"(\d+)\s*(?:of|/|out of)\s*(\d+)\s*(?:carriers|individuals|relatives|heterozygotes)", re.IGNORECASE),
    # "12 affected ... 8 unaffected"
    re.compile(r"(\d+)\s*affected\D{0,40}?(\d+)\s*unaffected", re.IGNORECASE),
]


def _mentions_any_variant(sent: str) -> bool:
    return bool(_PROT_RE.search(sent) or _CDNA_RE.search(sent))


def _parse_counts(sent: str):
    for pat in _COUNT_PATTERNS:
        m = pat.search(sent)
        if not m:
            continue
        a = float(m.group(1))
        second = float(m.group(2))
        if "unaffected" in pat.pattern:
            return a, second
        return a, max(second - a, 0.0)  # "a of n" -> n is the total
    return None


def regex_carrier_extractor(
    record: LiteratureRecord, tokens: Sequence[str]
) -> List[ExtractionResult]:
    """Dependency-free fallback extractor.

    Splits the document into sentences and tracks the "current variant context":
    the most recently mentioned variant. Count patterns are attributed to that
    context, so counts stated one sentence after the variant is named are still
    captured, while counts attached to a *different* variant are rejected. This
    is the variant-matching gate.
    """

    results: List[ExtractionResult] = []
    sentences = re.split(r"(?<=[.!?])\s+", record.text)
    context = None  # "target", "other", or None
    for sent in sentences:
        if _variant_mentioned(sent, tokens):
            context = "target"
        elif _mentions_any_variant(sent):
            context = "other"
        counts = _parse_counts(sent)
        if counts is None or context != "target":
            continue
        affected, unaffected = counts
        results.append(
            ExtractionResult(
                variant_token=next(iter(tokens)),
                affected=affected,
                unaffected=unaffected,
                phenotype=record.metadata.get("phenotype"),
                quote=sent.strip(),
            )
        )
    return results


class LiteratureCountAdapter(CountAdapter):
    name = "literature"

    def __init__(
        self,
        records: Sequence[LiteratureRecord],
        extractor: Extractor = regex_carrier_extractor,
        variant_table: Optional[pd.DataFrame] = None,
        phenotype_ontology: Optional[Dict[str, str]] = None,
    ):
        self._records = list(records)
        self._extractor = extractor
        if variant_table is None:
            variant_table = load_variant_labels()
        self._variant_table = variant_table.set_index("variant_id", drop=False)
        self._ontology = {k.lower(): v for k, v in (phenotype_ontology or {}).items()}

    def _tokens_for(self, variant_id: str) -> List[str]:
        hgvs = ""
        gene = ""
        if variant_id in self._variant_table.index:
            row = self._variant_table.loc[variant_id]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            hgvs = str(row.get("hgvs", ""))
            gene = str(row.get("gene", ""))
        tokens = normalize_variant_tokens(hgvs=hgvs, variant_id=variant_id)
        return tokens, gene

    def _normalize_phenotype(self, phenotype: Optional[str]) -> Optional[str]:
        if not phenotype:
            return None
        return self._ontology.get(phenotype.lower(), phenotype)

    def fetch(self, variant_id: str) -> Optional[CarrierCounts]:
        tokens, gene = self._tokens_for(variant_id)
        if not tokens:
            return None
        affected = 0.0
        unaffected = 0.0
        provenance: List[str] = []
        found = False
        for record in self._records:
            if gene and record.gene and record.gene.upper() != gene.upper():
                continue
            for res in self._extractor(record, tokens):
                found = True
                affected += res.affected
                unaffected += res.unaffected
                pheno = self._normalize_phenotype(res.phenotype)
                pheno_str = f" [{pheno}]" if pheno else ""
                provenance.append(f"{record.doc_id}{pheno_str}: \"{res.quote}\"")
        if not found:
            return None
        return CarrierCounts(
            variant_id=variant_id,
            affected=affected,
            unaffected=unaffected,
            source=self.name,
            provenance=provenance,
        )
