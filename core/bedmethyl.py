"""bedMethyl helpers: count validation and pileup mod-base expansion."""
from __future__ import annotations

import gzip
from pathlib import Path


def pileup_modified_bases(modified_bases: tuple[str, ...]) -> tuple[str, ...]:
    """Return ``modified_bases`` augmented for a faithful modkit pileup.

    modkit DMR requires every modification type present in the modBAM to appear
    in the pileup table. When the analysis targets 5mC only, we still tabulate
    5hmC so ``N_other`` counts are consistent; :func:`core.pileup.fold_hmc`
    folds confident 5hmC into canonical before DMR.
    """
    bases = list(modified_bases)
    if "5mC" in bases and "5hmC" not in bases:
        bases.append("5hmC")
    return tuple(bases)


def _parse_counts(line: str) -> tuple[int, int, int, int] | None:
    if line.startswith("#") or not line.strip():
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 14:
        return None
    try:
        valid = int(fields[9])
        n_mod = int(fields[11])
        n_canonical = int(fields[12])
        n_other = int(fields[13])
    except ValueError:
        return None
    return valid, n_mod, n_canonical, n_other


def validate_bedmethyl(
    path: str | Path,
    *,
    max_examples: int = 5,
) -> tuple[int, int, list[str]]:
    """Check that ``valid_coverage == N_mod + N_canonical + N_other`` per row.

    Returns ``(n_rows, n_invalid, example_lines)``. ``example_lines`` are raw
    bedMethyl rows (truncated) for the first invalid positions.
    """
    source = Path(path)
    opener = gzip.open if str(source).endswith(".gz") else open
    n_rows = n_invalid = 0
    examples: list[str] = []
    with opener(source, "rt") as handle:
        for line in handle:
            parsed = _parse_counts(line)
            if parsed is None:
                continue
            fields = line.rstrip("\n").split("\t")
            n_rows += 1
            valid, n_mod, n_canonical, n_other = parsed
            if valid != n_mod + n_canonical + n_other:
                n_invalid += 1
                if len(examples) < max_examples:
                    start = fields[1] if len(fields) > 1 else "?"
                    examples.append(
                        f"{fields[0]}:{start} valid={valid} "
                        f"mod={n_mod} canon={n_canonical} other={n_other}"
                    )
    return n_rows, n_invalid, examples


def bedmethyl_has_rows(path: str | Path) -> bool:
    """True when the bedMethyl file contains at least one data row."""
    source = Path(path)
    if not source.exists() or source.stat().st_size == 0:
        return False
    opener = gzip.open if str(source).endswith(".gz") else open
    with opener(source, "rt") as handle:
        for line in handle:
            if _parse_counts(line) is not None:
                return True
    return False
