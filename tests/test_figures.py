from pathlib import Path

import pandas as pd
import pysam

from core import figures


def _write_pileup(path: Path, rows: list[tuple[int, float, int]]) -> Path:
    plain = path.with_suffix(".plain.bed")
    with plain.open("w") as handle:
        for pos, pct, cov in rows:
            handle.write("\t".join([
                "chr3", str(pos), str(pos + 1), "m", str(cov), "+",
                str(pos), str(pos + 1), "0,0,0", str(cov), f"{pct}",
                "0", "0", "0", "0", "0", "0", "0",
            ]) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)
    plain.unlink()
    return path


def test_gene_locus_plot_writes_png(tmp_path: Path):
    positions = list(range(57192900, 57193100, 20))
    proband = _write_pileup(tmp_path / "P.bed.gz", [(p, 13.0, 25) for p in positions])
    r1 = _write_pileup(tmp_path / "M.bed.gz", [(p, 12.0, 22) for p in positions])
    r2 = _write_pileup(tmp_path / "B.bed.gz", [(p, 12.0, 20) for p in positions])

    out = figures.gene_locus_plot(
        "HESX1", "chr3", 57192837, 57232606, (57192837, 57194837),
        {"Proband": proband, "Mother": r1, "Brother": r2},
        tmp_path / "gene_HESX1.png",
        min_cov=10, subtitle="promoter P 13% M 12% B 12%",
        colors={"Proband": figures.ROLE_COLORS[0], "Mother": figures.ROLE_COLORS[1],
                "Brother": figures.ROLE_COLORS[2]},
    )
    assert out is not None and out.exists() and out.stat().st_size > 0


def test_gene_locus_plot_returns_none_without_data(tmp_path: Path):
    empty = _write_pileup(tmp_path / "P.bed.gz", [])
    out = figures.gene_locus_plot(
        "HESX1", "chr9", 1, 100, (1, 50), {"Proband": empty},
        tmp_path / "none.png",
    )
    assert out is None


def test_targeted_heatmap(tmp_path: Path):
    scored = pd.DataFrame([
        {"gene": "HESX1", "region": "promoter", "candidate": True,
         "proband_meth": 13.0, "Mother_meth": 12.0, "Brother_meth": 12.0},
        {"gene": "SOX2", "region": "promoter", "candidate": False,
         "proband_meth": 80.0, "Mother_meth": 78.0, "Brother_meth": 79.0},
    ])
    out = figures.targeted_heatmap(
        scored, "promoter", ["proband_meth", "Mother_meth", "Brother_meth"],
        ["Proband", "Mother", "Brother"], tmp_path / "heatmap.png",
    )
    assert out is not None and out.exists()


def test_karyotype_and_histogram(tmp_path: Path):
    candidates = pd.DataFrame([
        {"chrom": "chr3", "start": 57192837, "end": 57232606, "mean_effect": 0.4,
         "imprinted_control": ""},
        {"chrom": "chr11", "start": 1999000, "end": 2003000, "mean_effect": -0.3,
         "imprinted_control": "H19/ICR1"},
    ])
    null = pd.DataFrame({"effect": [0.02, -0.05, 0.1, 0.03, -0.2]})
    karyo = figures.karyotype_plot(candidates, tmp_path / "karyo.png")
    hist = figures.effect_histogram(candidates, null, tmp_path / "hist.png")
    assert karyo.exists() and hist.exists()


def test_karyotype_handles_empty(tmp_path: Path):
    out = figures.karyotype_plot(pd.DataFrame(columns=["chrom", "start", "end", "mean_effect"]),
                                 tmp_path / "empty.png")
    assert out.exists()


def test_karyotype_uses_gene_labels(tmp_path: Path):
    candidates = pd.DataFrame([
        {"chrom": "chr3", "start": 100, "end": 200, "mean_effect": 0.4,
         "promoter_of": "HESX1", "genes": "HESX1,OTHER", "imprinted_control": ""},
    ])
    out = figures.karyotype_plot(candidates, tmp_path / "k.png")
    assert out.exists()


def test_sweep_plot(tmp_path: Path):
    thresholds = [0.1, 0.2, 0.3, 0.4]
    series = [
        {"label": "P (proband)", "counts": [19, 16, 16, 10], "color": "#C0392B", "linestyle": "-"},
        {"label": "B (sibling control)", "counts": [11, 11, 10, 6], "color": "#2471A3", "linestyle": "-"},
        {"label": "M (parent)", "counts": [7, 7, 7, 6], "color": "#F39C12", "linestyle": "--"},
    ]
    out = figures.sweep_plot(thresholds, series, tmp_path / "sweep.png",
                             ylabel="private DMRs (autosomes)", title="sweep")
    assert out.exists() and out.stat().st_size > 0
