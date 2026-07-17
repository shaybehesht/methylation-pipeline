"""Sliding-window methylation profiles for the interactive (zoomable) viewer."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from core.targeted import load_region


def region_profile(
    pileups: dict[str, Path], chrom: str, start: int, end: int,
    *, window: int = 20, min_cov: int = 1,
) -> pd.DataFrame:
    """Return a tidy per-CpG frame with a sliding-window mean per sample.

    Columns: ``sample, pos, pct, smooth``. ``window`` is the rolling window in
    CpG sites (centered), matching the reference "sliding window mean" profile.
    """
    frames = []
    for label, path in pileups.items():
        methylation = load_region(path, chrom, int(start), int(end), min_cov)
        if not methylation:
            continue
        positions = sorted(methylation)
        frame = pd.DataFrame({
            "sample": label,
            "pos": positions,
            "pct": [methylation[pos] for pos in positions],
        })
        span = max(1, int(window))
        frame["smooth"] = frame["pct"].rolling(span, center=True, min_periods=1).mean()
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["sample", "pos", "pct", "smooth"])
    return pd.concat(frames, ignore_index=True)


def parse_region(text: str) -> tuple[str, int, int] | None:
    """Parse ``chrom:start-end`` (commas allowed) into ``(chrom, start, end)``."""
    try:
        chrom, span = text.strip().split(":")
        low, high = span.replace(",", "").split("-")
        start, end = int(low), int(high)
    except (ValueError, AttributeError):
        return None
    if not chrom or end <= start:
        return None
    return chrom, start, end
