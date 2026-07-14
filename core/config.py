"""Validated run configuration and trio comparison derivation."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any


class Sex(str, Enum):
    FEMALE = "F"
    MALE = "M"


class Role(str, Enum):
    PROBAND = "proband"
    RELATIVE = "relative"


class Affection(str, Enum):
    AFFECTED = "affected"
    UNAFFECTED = "unaffected"
    UNKNOWN = "unknown"


class Relationship(str, Enum):
    MOTHER = "mother"
    FATHER = "father"
    SIBLING = "sibling"
    OTHER = "other"


@dataclass(frozen=True)
class Sample:
    label: str
    bam_path: str
    sex: Sex
    role: Role
    relationship: Relationship | None = None
    affection: Affection | None = None

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("Sample labels cannot be empty")
        if not self.bam_path.strip():
            raise ValueError(f"BAM path is required for {self.label}")


@dataclass(frozen=True)
class Comparison:
    name: str
    left: Sample
    right: Sample
    kind: str

    def valid_chromosome(self, chromosome: str) -> bool:
        chrom = chromosome.lower().removeprefix("chr")
        if chrom in {"m", "mt"}:
            return False
        if chrom == "x":
            return self.left.sex == self.right.sex == Sex.FEMALE
        if chrom == "y":
            return self.left.sex == self.right.sex == Sex.MALE
        return True


@dataclass
class RegionConfig:
    mode: str = "chromosomes"
    chromosomes: list[str] = field(default_factory=lambda: ["chr1", "chr2", "chr11", "chr15"])
    genes: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.mode not in {"whole_genome", "chromosomes", "targeted"}:
            raise ValueError(f"Unsupported region mode: {self.mode}")
        if self.mode == "chromosomes" and not self.chromosomes:
            raise ValueError("Select at least one chromosome")
        if self.mode == "targeted" and not self.genes:
            raise ValueError("Provide at least one target gene")


@dataclass
class TrioConfig:
    samples: list[Sample]
    reference_fasta: str
    output_dir: str = "runs/latest"
    regions: RegionConfig = field(default_factory=RegionConfig)
    thresholds: dict[str, float | int] = field(default_factory=dict)
    phased_vcf: str = ""
    assembly: str = ""
    modified_bases: list[str] = field(default_factory=lambda: ["5mC"])
    combine_strands: bool = True

    def __post_init__(self) -> None:
        if len(self.samples) != 3:
            raise ValueError("Exactly three samples are required")
        labels = [sample.label for sample in self.samples]
        if len(set(labels)) != 3:
            raise ValueError("Sample labels must be unique")
        if sum(sample.role == Role.PROBAND for sample in self.samples) != 1:
            raise ValueError("Exactly one sample must be the proband")
        if not self.reference_fasta.strip():
            raise ValueError("Reference FASTA is required")

    @property
    def proband(self) -> Sample:
        return next(sample for sample in self.samples if sample.role == Role.PROBAND)

    @property
    def relatives(self) -> list[Sample]:
        return [sample for sample in self.samples if sample.role == Role.RELATIVE]

    def comparisons(self) -> list[Comparison]:
        p, (r1, r2) = self.proband, self.relatives
        return [
            Comparison(f"{p.label}_vs_{r1.label}", p, r1, "proband"),
            Comparison(f"{p.label}_vs_{r2.label}", p, r2, "proband"),
            Comparison(f"{r1.label}_vs_{r2.label}", r1, r2, "null"),
        ]

    def analysis_design(self) -> str:
        statuses = [sample.affection for sample in self.relatives]
        if statuses.count(Affection.UNAFFECTED) == 2:
            return "proband_specific"
        if statuses.count(Affection.AFFECTED) == 1 and statuses.count(Affection.UNAFFECTED) == 1:
            return "phenotype_segregation"
        if statuses.count(Affection.AFFECTED) == 2:
            return "no_unaffected_control"
        return "phenotype_unknown"

    def evidence_status(self) -> dict[str, str]:
        relationships = {sample.relationship for sample in self.relatives}
        both_parents = {Relationship.MOTHER, Relationship.FATHER} <= relationships
        return {
            "phenotype": self.analysis_design(),
            "parent_of_origin": (
                "inputs_available" if both_parents and self.phased_vcf
                else "needs_both_parents_and_phased_vcf"
            ),
            "mqtl": "phased_vcf_available" if self.phased_vcf else "needs_phased_vcf",
        }

    def caveats(self) -> list[str]:
        sexes = {sample.sex for sample in self.relatives}
        notes = ["Three samples are a family comparison, not a population cohort."]
        design = self.analysis_design()
        if design == "phenotype_segregation":
            notes.append(
                "One affected and one unaffected relative were provided; candidates are "
                "ranked for similarity within affected samples and difference from the unaffected sample."
            )
        elif design == "phenotype_unknown":
            notes.append(
                "Relative clinical status is incomplete; phenotype segregation cannot be assessed."
            )
        elif design == "no_unaffected_control":
            notes.append(
                "All relatives are affected; no unaffected family control is available."
            )
        evidence = self.evidence_status()
        if not self.phased_vcf:
            notes.append("No phased VCF was provided; mQTL and parent-of-origin effects are not resolved.")
        elif evidence["parent_of_origin"] != "inputs_available":
            notes.append("A phased VCF is present, but both identified parents are needed to assign parental origin.")
        if Sex.MALE not in sexes:
            notes.append(
                "No male relative is present; paternal-allele mQTL effects cannot be "
                "separated and are an important false-positive mode."
            )
        if Sex.FEMALE not in sexes:
            notes.append(
                "No female relative is present; maternal-allele mQTL effects cannot be "
                "separated and are an important false-positive mode."
            )
        return notes

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def ensure_output_dir(self) -> Path:
        path = Path(self.output_dir).expanduser().resolve()
        path.mkdir(parents=True, exist_ok=True)
        return path
