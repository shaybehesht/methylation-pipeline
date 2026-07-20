"""Maximum credible allele frequency bound (Whiffin/Ware, Genet Med 2017).

For a dominant, rare-allele disease the carrier frequency is approximately
``2 * AF``. Bounding penetrance from above:

    p <= K * f_allelic * f_genetic / (2 * AF)

where ``K`` is disease prevalence, ``f_allelic`` the maximum fraction of cases
attributable to a single allele (allelic heterogeneity), and ``f_genetic`` the
maximum fraction attributable to the gene (genetic heterogeneity). Inverting it
gives the maximum credible allele frequency a genuinely pathogenic, penetrant
variant could have:

    AF_max = K * f_allelic * f_genetic / (2 * p)

An allele observed *above* ``AF_max`` in the population is "too common" to be
highly penetrant - the exact signature the gene model exploits in its
allele-frequency-spectrum feature. Reused here as a bound, not claimed as novel.
"""

from __future__ import annotations


def max_credible_af(
    prevalence: float,
    penetrance: float,
    f_allelic: float = 1.0,
    f_genetic: float = 1.0,
    inheritance: str = "dominant",
) -> float:
    """Maximum credible allele frequency for a pathogenic, penetrant variant."""

    _validate(prevalence, penetrance, f_allelic, f_genetic)
    max_carrier_freq = prevalence * f_allelic * f_genetic / penetrance
    if inheritance == "dominant":
        return float(max_carrier_freq / 2.0)
    if inheritance == "recessive":
        # Genotype frequency ~ AF^2 for a homozygous recessive model.
        return float(max_carrier_freq ** 0.5)
    raise ValueError("inheritance must be 'dominant' or 'recessive'")


def penetrance_upper_bound_from_af(
    allele_frequency: float,
    prevalence: float,
    f_allelic: float = 1.0,
    f_genetic: float = 1.0,
    inheritance: str = "dominant",
) -> float:
    """Upper bound on penetrance implied by an observed population allele frequency.

    Returns a value in ``[0, 1]`` (clamped): if the AF is small the bound is
    uninformative (1.0); if the allele is common the bound drops below 1.
    """

    if allele_frequency < 0:
        raise ValueError("allele_frequency must be non-negative")
    if allele_frequency == 0:
        return 1.0
    _validate(prevalence, 1.0, f_allelic, f_genetic)
    if inheritance == "dominant":
        carrier_freq = 2.0 * allele_frequency
    elif inheritance == "recessive":
        carrier_freq = allele_frequency ** 2
    else:
        raise ValueError("inheritance must be 'dominant' or 'recessive'")
    bound = prevalence * f_allelic * f_genetic / carrier_freq
    return float(min(1.0, max(0.0, bound)))


def _validate(prevalence: float, penetrance: float, f_allelic: float, f_genetic: float) -> None:
    if not 0 < prevalence <= 1:
        raise ValueError("prevalence must be in (0, 1]")
    if not 0 < penetrance <= 1:
        raise ValueError("penetrance must be in (0, 1]")
    if not 0 < f_allelic <= 1:
        raise ValueError("f_allelic must be in (0, 1]")
    if not 0 < f_genetic <= 1:
        raise ValueError("f_genetic must be in (0, 1]")
