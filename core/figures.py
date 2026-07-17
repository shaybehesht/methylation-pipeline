"""Publication-style figures matching the reference scripts (03/07/11).

* ``gene_locus_plot`` — per-gene % CpG methylation across a region for every
  sample, with the promoter highlighted and a coverage track (script 07).
* ``targeted_heatmap`` — promoter/body methylation heatmap with candidates boxed
  (script 07).
* ``karyotype_plot`` / ``effect_histogram`` — genome-wide proband-private DMRs
  and the family null (scripts 03/11).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from core.targeted import load_region  # noqa: E402

# Proband red, relatives orange/blue — matches the reference figures.
ROLE_COLORS = ["#C0392B", "#F39C12", "#2471A3"]
HYPER, HYPO, NULLC = "#C0392B", "#2471A3", "#BDC3C7"

CHROM_LEN = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895,
}
KARYO_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def _coverage(pileup_gz: str | Path, chrom: str, start: int, end: int):
    import pysam

    positions, coverage = [], []
    try:
        with pysam.TabixFile(str(pileup_gz)) as tabix:
            for line in tabix.fetch(chrom, max(0, start), end):
                fields = line.split("\t")
                if len(fields) < 10:
                    continue
                positions.append(int(fields[1]))
                coverage.append(int(fields[9]))
    except (OSError, ValueError):
        pass
    return positions, coverage


def gene_locus_plot(
    gene: str, chrom: str, start: int, end: int, promoter: tuple[int, int],
    pileups: dict[str, Path], out_png: str | Path, *, min_cov: int = 10,
    subtitle: str = "", colors: dict[str, str] | None = None,
) -> Path | None:
    """Per-gene methylation + coverage plot; returns the path, or None if no data."""
    colors = colors or {}
    fig, (ax, cover_ax) = plt.subplots(
        2, 1, figsize=(13, 6.2), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    drew = False
    for label, path in pileups.items():
        methylation = load_region(path, chrom, start, end, min_cov=1)
        color = colors.get(label, "#333333")
        if methylation:
            drew = True
            positions = np.array(sorted(methylation))
            percents = np.array([methylation[pos] for pos in positions])
            ax.plot(positions, percents, "o", ms=2.4, color=color, alpha=0.45)
            window = max(3, len(positions) // 25)
            smooth = pd.Series(percents).rolling(window, center=True, min_periods=1).mean()
            ax.plot(positions, smooth, lw=2.4, color=color, label=label)
        cov_pos, cov_val = _coverage(path, chrom, start, end)
        if cov_pos:
            cover_ax.plot(cov_pos, cov_val, lw=0.9, color=color, alpha=0.7)
    if not drew:
        plt.close(fig)
        return None

    tss_lo, tss_hi = promoter
    ax.axvspan(tss_lo, tss_hi, color="#F9E79F", alpha=0.5, zorder=0, label="promoter")
    ax.set_ylim(-3, 103)
    ax.set_ylabel("% CpG methylation")
    ax.legend(frameon=False, ncol=4, fontsize=9, loc="upper left")
    ax.set_title(f"{gene}   {chrom}:{start:,}-{end:,}   {subtitle}", fontsize=12)
    cover_ax.axvspan(tss_lo, tss_hi, color="#F9E79F", alpha=0.5, zorder=0)
    cover_ax.axhline(min_cov, color="#95A5A6", ls="--", lw=0.8)
    cover_ax.set_ylabel("coverage")
    cover_ax.set_xlabel(f"{chrom} (bp)")
    for axis in (ax, cover_ax):
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=190)
    plt.close(fig)
    return Path(out_png)


def targeted_heatmap(
    scored: pd.DataFrame, kind: str, meth_columns: list[str], sample_labels: list[str],
    out_png: str | Path,
) -> Path | None:
    """Heatmap of promoter/body methylation across samples; candidates boxed."""
    subset = scored[(scored["region"] == kind) & scored[meth_columns[0]].notna()]
    if subset.empty:
        return None
    subset = subset.sort_values("gene")
    matrix = subset[meth_columns].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(6.4, max(6, 0.23 * len(subset))))
    image = ax.imshow(matrix, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100)
    ax.set_xticks(range(len(sample_labels)))
    ax.set_xticklabels(sample_labels, fontsize=9)
    ax.set_yticks(range(len(subset)))
    ax.set_yticklabels(list(subset["gene"]), fontsize=7)
    if "candidate" in subset.columns:
        for row, is_candidate in enumerate(subset["candidate"]):
            if is_candidate:
                ax.add_patch(plt.Rectangle(
                    (-0.5, row - 0.5), len(sample_labels), 1, fill=False,
                    edgecolor="k", lw=2.2,
                ))
    ax.set_title(f"{kind} methylation\n(black box = proband-specific)", fontsize=11)
    fig.colorbar(image, ax=ax, label="% CpG methylation", fraction=0.04)
    fig.tight_layout()
    fig.savefig(out_png, dpi=190)
    plt.close(fig)
    return Path(out_png)


def sweep_plot(thresholds, series: list[dict], out_png: str | Path, *, ylabel: str, title: str) -> Path:
    """Private-DMR counts vs |effect| threshold for the proband and relatives."""
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    for entry in series:
        ax.plot(
            thresholds, entry["counts"], marker="o", ms=4,
            lw=entry.get("linewidth", 2.4), color=entry["color"],
            linestyle=entry.get("linestyle", "-"), alpha=entry.get("alpha", 1.0),
            label=entry["label"],
        )
    ax.set_xlabel("|effect size| threshold")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return Path(out_png)


def karyotype_plot(
    candidates: pd.DataFrame, out_png: str | Path, *, title: str = "Proband-specific DMRs",
    effect_column: str = "mean_effect",
    label_columns: tuple[str, ...] = ("promoter_of", "genes", "imprinted_control"),
) -> Path:
    """Karyotype-style lollipop plot of proband-private DMRs (scripts 03/11)."""
    from matplotlib.patches import Rectangle

    chroms = [c for c in KARYO_ORDER if c in CHROM_LEN]
    fig, ax = plt.subplots(figsize=(15, 9))
    maxlen = max(CHROM_LEN[c] for c in chroms)
    for index, chrom in enumerate(chroms):
        y = len(chroms) - index
        ax.add_patch(Rectangle(
            (0, y - 0.11), CHROM_LEN[chrom], 0.22,
            facecolor="#ECF0F1", edgecolor="#95A5A6", lw=0.6, zorder=1,
        ))
        if candidates.empty or effect_column not in candidates.columns:
            continue
        for _, region in candidates[candidates["chrom"] == chrom].iterrows():
            effect = float(region[effect_column])
            height = min(abs(effect), 1.0) * 0.85
            sign = 1 if effect > 0 else -1
            x = (int(region["start"]) + int(region["end"])) / 2
            color = HYPER if effect > 0 else HYPO
            ax.plot([x, x], [y, y + sign * height], color=color, lw=1.2, alpha=0.9, zorder=3)
            ax.plot([x], [y + sign * height], "o", ms=3.4, color=color, zorder=4)
            label = ""
            for column in label_columns:
                value = str(region.get(column, "") or "")
                if value:
                    label = value.split(",")[0]
                    break
            if label:
                ax.annotate(
                    label, (x, y + sign * height), xytext=(0, 7 * sign),
                    textcoords="offset points", ha="center", fontsize=6, color="#7D3C98",
                )
    ax.set_yticks(range(1, len(chroms) + 1))
    ax.set_yticklabels(list(reversed(chroms)), fontsize=8)
    ax.set_xlim(-maxlen * 0.01, maxlen * 1.02)
    ax.set_ylim(0.2, len(chroms) + 1.1)
    ax.set_xlabel("position (Mb)")
    ax.set_xticks(np.arange(0, maxlen, 25_000_000))
    ax.set_xticklabels([str(int(v / 1e6)) for v in np.arange(0, maxlen, 25_000_000)], fontsize=8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.plot([], [], color=HYPER, lw=2, label="hypermethylated in proband")
    ax.plot([], [], color=HYPO, lw=2, label="hypomethylated in proband")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title(f"{title}  (n={0 if candidates.empty else len(candidates)})", fontsize=13, pad=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return Path(out_png)


def effect_histogram(
    candidates: pd.DataFrame, null_segments: pd.DataFrame, out_png: str | Path,
) -> Path:
    """Proband hit effect sizes (lines) over the relative–relative null (grey)."""
    fig, ax = plt.subplots(figsize=(8, 5.2))
    if not null_segments.empty and "effect" in null_segments.columns:
        values = null_segments["effect"].dropna()
        if len(values):
            ax.hist(values, bins=40, color=NULLC, alpha=0.75, density=True,
                    label="relative vs relative (null)")
    if not candidates.empty and "mean_effect" in candidates.columns:
        for effect in candidates["mean_effect"]:
            ax.axvline(effect, color=HYPER if effect > 0 else HYPO, lw=1.0, alpha=0.8)
    ax.set_xlabel("effect size")
    ax.set_ylabel("density")
    ax.set_title("Proband-private DMRs (lines) vs the family null (grey)", fontsize=12)
    ax.legend(frameon=False, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return Path(out_png)
