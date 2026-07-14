from pathlib import Path

import pandas as pd
import pytest

from core.dmr import build_command as dmr_command
from core.pileup import build_command as pileup_command, normalize_counts
from core.reasoning import write_html_report
from core.thresholds import REGISTRY, defaults, validate, widget_values


def test_modkit_commands_are_reproducible():
    pileup = pileup_command("a.bam", "a.bed.gz", "ref.fa", include_bed="scope.bed")
    assert pileup[:3] == ["modkit", "pileup", "a.bam"]
    assert "--bgzf" in pileup and "--include-bed" in pileup
    # modkit >=0.6 requires --modified-bases when --cpg is used.
    assert "--modified-bases" in pileup
    assert pileup[pileup.index("--modified-bases") + 1] == "5mC"
    assert pileup.index("--modified-bases") < pileup.index("--cpg")
    dmr = dmr_command("a.bed.gz", "b.bed.gz", "out.tsv", "ref.fa", regions="scope.bed")
    assert "--regions-bed" in dmr and "--header" in dmr and "--base" in dmr


def test_pileup_supports_multiple_modified_bases():
    pileup = pileup_command("a.bam", "a.bed.gz", "ref.fa", modified_bases=("5mC", "5hmC"))
    index = pileup.index("--modified-bases")
    assert pileup[index + 1:index + 3] == ["5mC", "5hmC"]
    assert pileup[index + 3] == "--cpg"


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
