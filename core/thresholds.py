"""Single source of truth for controls and their explanations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class Threshold:
    label: str
    default: float | int
    widget: Literal["slider", "number"]
    minimum: float | int
    maximum: float | int
    step: float | int
    rationale: str


REGISTRY: dict[str, Threshold] = {
    "null_percentile": Threshold(
        "Null effect percentile", 99, "slider", 90, 100, 0.5,
        "Uses the extreme tail of relative-versus-relative effects as the family's empirical noise floor.",
    ),
    "min_sites": Threshold(
        "Minimum CpG sites", 5, "number", 1, 100, 1,
        "Multiple CpGs reduce the chance that one noisy site defines a region.",
    ),
    "max_pval": Threshold(
        "Maximum p-value", 0.01, "number", 0.0001, 0.1, 0.001,
        "A stringent unadjusted screen limits candidates in this n=3 exploratory design.",
    ),
    "min_valid_coverage": Threshold(
        "Minimum valid coverage", 10, "number", 1, 100, 1,
        "Ten informative reads is a practical floor for stable methylation fractions.",
    ),
    "filter_threshold": Threshold(
        "Modkit probability threshold", 0.7, "slider", 0.5, 1.0, 0.01,
        "Calls below 0.7 are treated as ambiguous rather than methylated or unmethylated.",
    ),
    "targeted_min_delta": Threshold(
        "Targeted minimum delta (percentage points)", 10, "number", 0, 100, 1,
        "Ten percentage points prioritizes changes large enough to inspect biologically.",
    ),
    "alpha": Threshold(
        "Targeted alpha", 0.05, "number", 0.001, 0.2, 0.001,
        "Conventional exploratory significance threshold for targeted Wilcoxon tests.",
    ),
    "promoter_pad": Threshold(
        "Promoter padding (bp)", 2000, "number", 0, 10000, 100,
        "Two kilobases captures common proximal promoter regulatory sequence.",
    ),
    "body_pad": Threshold(
        "Gene-body padding (bp)", 5000, "number", 0, 50000, 500,
        "Five kilobases includes nearby regulatory context while keeping targeted runs compact.",
    ),
}


def defaults() -> dict[str, float | int]:
    return {key: value.default for key, value in REGISTRY.items()}


def widget_values(
    spec: Threshold, current: float | int
) -> tuple[float | int, float | int, float | int, float | int]:
    """Return bounds, value, and step with Streamlit-compatible matching types."""
    values = (spec.minimum, spec.maximum, current, spec.step)
    numeric_type = float if any(isinstance(value, float) for value in values) else int
    return tuple(numeric_type(value) for value in values)


def validate(values: dict[str, float | int]) -> dict[str, float | int]:
    result = defaults()
    for key, value in values.items():
        if key not in REGISTRY:
            raise ValueError(f"Unknown threshold: {key}")
        spec = REGISTRY[key]
        if not spec.minimum <= value <= spec.maximum:
            raise ValueError(f"{key} must be between {spec.minimum} and {spec.maximum}")
        result[key] = value
    return result
