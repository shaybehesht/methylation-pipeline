from pathlib import Path

import pandas as pd
import pytest

from core.dmr import build_command as dmr_command
from core.pileup import build_command as pileup_command, normalize_counts
from core.reasoning import write_html_report
from core.thresholds import REGISTRY, defaults, validate, widget_values


def test_modkit_commands_are_reproducible():
    pileup = pileup_command("a.bam", "a.bed.gz", "ref.fa", region="chr1")
    assert pileup[:3] == ["modkit", "pileup", "a.bam"]
    assert "--cpg" in pileup and "--combine-strands" in pileup
    assert "--region" in pileup and pileup[pileup.index("--region") + 1] == "chr1"
    # modkit >=0.6 requires --modified-bases when --cpg is used.
    assert "--modified-bases" in pileup
    assert pileup[pileup.index("--modified-bases") + 1] == "5mC"
    assert pileup.index("--cpg") < pileup.index("--modified-bases")
    # We compress/index ourselves, so --bgzf is not passed to modkit.
    assert "--bgzf" not in pileup
    # --chunk-size is not a real modkit pileup flag (it rejects it with exit 2);
    # --interval-size is the supported memory control.
    assert "--chunk-size" not in pileup
    dmr = dmr_command("a.bed.gz", "b.bed.gz", "out.tsv", "ref.fa", segment="seg.bed")
    assert "--segment" in dmr and "--header" in dmr and "--base" in dmr and "--force" in dmr


def test_pileup_can_disable_combine_strands_and_take_multiple_mods():
    pileup = pileup_command(
        "a.bam", "a.bed.gz", "ref.fa", combine_strands=False, modified_bases=("5mC", "5hmC")
    )
    assert "--combine-strands" not in pileup
    index = pileup.index("--modified-bases")
    assert pileup[index + 1:index + 3] == ["5mC", "5hmC"]
    assert pileup[index + 3] == "--filter-threshold"


def test_general_worker_path_uses_motif_without_banned_flags():
    # PacBio HiFi route: general workers (--motif CG 0) instead of the optimized
    # --cpg workers, and no --modified-bases. Never emit flags that bioconda
    # modkit 0.6.x rejects (--chunk-size, --force-allow-implicit).
    pileup = pileup_command(
        "a.bam", "a.bed.gz", "ref.fa", region="chr14",
        combine_strands=False, use_general_workers=True,
    )
    assert "--cpg" not in pileup
    assert "--modified-bases" not in pileup
    assert "--force-allow-implicit" not in pileup
    assert "--chunk-size" not in pileup
    motif = pileup.index("--motif")
    assert pileup[motif + 1:motif + 3] == ["CG", "0"]
    assert pileup[pileup.index("--region") + 1] == "chr14"


def test_pileup_commands_never_emit_banned_flags():
    for kwargs in (
        {},
        {"combine_strands": False, "use_general_workers": True},
        {"combine_strands": False, "modified_bases": ("5mC", "5hmC")},
    ):
        command = pileup_command(
            "a.bam", "a.bed.gz", "ref.fa", region="chr14",
            threads=8, interval_size=100_000, **kwargs,
        )
        for banned in ("--chunk-size", "--force-allow-implicit"):
            assert banned not in command
        # Every flag we emit must look like a real option (starts with --)
        # or a positional / value. Spot-check the known-good options.
        for flag in ("--ref", "--region", "--filter-threshold", "--threads",
                     "--interval-size", "--suppress-progress"):
            assert flag in command


def test_n_other_is_recomputed():
    source = pd.DataFrame({
        "N_valid_cov": [10, 5], "N_mod": [6, 1],
        "N_canonical": [3, 2], "N_other": [99, 99],
    })
    fixed = normalize_counts(source)
    assert fixed["N_other"].tolist() == [1, 2]


def test_threshold_validation():
    assert defaults()["null_percentile"] == 99
    with pytest.raises(ValueError, match="between"):
        validate({"alpha": 2})


def test_widget_values_have_matching_numeric_types():
    for spec in REGISTRY.values():
        values = widget_values(spec, spec.default)
        assert len({type(value) for value in values}) == 1
    assert widget_values(REGISTRY["null_percentile"], 99) == (90.0, 100.0, 99.0, 0.5)


def test_report_is_self_contained(tmp_path: Path):
    figure = tmp_path / "plot.png"
    figure.write_bytes(b"\x89PNG\r\n\x1a\n")
    output = write_html_report(
        tmp_path / "report.html", "Test", {"verdict": "MARGINAL"},
        "Reasoning", pd.DataFrame([{"rank": 1}]), figure,
    )
    text = output.read_text()
    assert "data:image/png;base64" in text
    assert "MARGINAL" in text
