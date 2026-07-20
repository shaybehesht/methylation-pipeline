"""Case/control + prevalence penetrance estimator (ADpenetrance-style).

For a disease with population prevalence ``K``:

    p = P(V|D) * K / [ P(V|D) * K + P(V|not D) * (1 - K) ]

where ``P(V|D)`` is the frequency of the variant genotype in cases and
``P(V|not D)`` its frequency in unaffected controls. This is Bayes' theorem
applied to ``P(affected | carrier)`` and is the core of the ADpenetrance model
(Wright/KCL). Reused here, not claimed as novel.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


def case_control_penetrance(
    freq_in_cases: float,
    freq_in_controls: float,
    prevalence: float,
) -> float:
    """Penetrance from variant frequency in cases vs controls and prevalence.

    Parameters
    ----------
    freq_in_cases:
        ``P(V | D)`` - carrier frequency among affected individuals.
    freq_in_controls:
        ``P(V | not D)`` - carrier frequency among unaffected individuals.
    prevalence:
        ``K`` - population prevalence of the disease.
    """

    for name, value in (
        ("freq_in_cases", freq_in_cases),
        ("freq_in_controls", freq_in_controls),
        ("prevalence", prevalence),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be in [0, 1], got {value}")
    numerator = freq_in_cases * prevalence
    denominator = numerator + freq_in_controls * (1.0 - prevalence)
    if denominator == 0:
        return 0.0
    return float(numerator / denominator)


@dataclass(frozen=True)
class CaseControlEstimate:
    penetrance: float
    ci_low: float
    ci_high: float


def case_control_penetrance_from_counts(
    carriers_in_cases: int,
    n_cases: int,
    carriers_in_controls: int,
    n_controls: int,
    prevalence: float,
    n_boot: int = 20000,
    level: float = 0.95,
    random_state: Optional[int] = 0,
) -> CaseControlEstimate:
    """Case/control penetrance with a Monte-Carlo credible interval.

    The two carrier frequencies are given Jeffreys ``Beta(0.5, 0.5)`` posteriors
    from the observed counts; we propagate their uncertainty through the Bayes
    formula by sampling.
    """

    if n_cases <= 0 or n_controls <= 0:
        raise ValueError("n_cases and n_controls must be positive")
    rng = np.random.default_rng(random_state)
    case_draws = rng.beta(carriers_in_cases + 0.5, n_cases - carriers_in_cases + 0.5, n_boot)
    ctrl_draws = rng.beta(
        carriers_in_controls + 0.5, n_controls - carriers_in_controls + 0.5, n_boot
    )
    num = case_draws * prevalence
    den = num + ctrl_draws * (1.0 - prevalence)
    samples = np.where(den > 0, num / den, 0.0)
    tail = (1.0 - level) / 2.0
    point = case_control_penetrance(
        carriers_in_cases / n_cases, carriers_in_controls / n_controls, prevalence
    )
    lo, hi = np.quantile(samples, [tail, 1.0 - tail])
    return CaseControlEstimate(penetrance=point, ci_low=float(lo), ci_high=float(hi))
