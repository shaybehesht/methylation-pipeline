"""Static figures suitable for both Streamlit and HTML reports."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def effect_plot(frame: pd.DataFrame, path: str | Path) -> Path:
    output = Path(path)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    if frame.empty:
        ax.text(0.5, 0.5, "No qualifying DMRs", ha="center", va="center")
        ax.set_axis_off()
    else:
        colors = ["#b2182b" if value > 0 else "#2166ac" for value in frame["effect_1"]]
        ax.scatter(frame["rank"], frame["mean_abs_effect"], c=colors, alpha=0.8)
        ax.set(xlabel="Candidate rank", ylabel="Mean absolute effect", title="Proband-specific DMR effects")
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output
