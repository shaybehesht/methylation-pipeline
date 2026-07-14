import json
from pathlib import Path

import pytest

from streamlit.testing.v1 import AppTest

SETUP_PAGE = str(Path(__file__).resolve().parents[1] / "app" / "pages" / "1_Setup.py")
RESULTS_PAGE = str(Path(__file__).resolve().parents[1] / "app" / "pages" / "5_Results.py")


def test_setup_page_renders_with_local_root(monkeypatch, tmp_path: Path):
    data_root = tmp_path / "data"
    (data_root / "cohort").mkdir(parents=True)
    (data_root / "cohort" / "proband.bam").write_bytes(b"")
    monkeypatch.setenv("METHYL_TRIO_DATA_ROOT", str(data_root))
    monkeypatch.setenv("METHYL_TRIO_REFERENCE_CACHE", str(tmp_path / "cache"))

    app = AppTest.from_file(SETUP_PAGE).run(timeout=30)
    assert not app.exception

    assembly_labels = [option for select in app.selectbox for option in select.options]
    assert any("hg38" in label or "GRCh38" in label for label in assembly_labels)
    assert any("Validate" in button.label for button in app.button)


def test_results_page_without_run_shows_hint(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("METHYL_TRIO_DATA_ROOT", str(tmp_path))
    monkeypatch.setenv("METHYL_TRIO_REFERENCE_CACHE", str(tmp_path / "cache"))
    app = AppTest.from_file(RESULTS_PAGE).run(timeout=30)
    assert not app.exception
    assert any("Run an analysis" in info.value for info in app.info)


def test_results_page_builds_complete_zip(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("METHYL_TRIO_DATA_ROOT", str(tmp_path))
    monkeypatch.setenv("METHYL_TRIO_REFERENCE_CACHE", str(tmp_path / "cache"))

    output = tmp_path / "run-1"
    output.mkdir()
    report = output / "report.html"
    report.write_text("<html></html>", encoding="utf-8")
    (output / "proband_specific_DMRs.tsv").write_text("rank\n1\n", encoding="utf-8")
    (output / "pipeline.log").write_text("done\n", encoding="utf-8")
    result = {
        "verdict": "MARGINAL",
        "candidate_count": 1,
        "ratio": 0.5,
        "reasoning": "Exploratory screen.",
        "evidence_status": {
            "phenotype": "phenotype_unknown",
            "parent_of_origin": "needs_both_parents_and_phased_vcf",
            "mqtl": "needs_phased_vcf",
        },
        "output": str(output),
        "report": str(report),
    }
    (output / "summary.json").write_text(json.dumps(result), encoding="utf-8")

    app = AppTest.from_file(RESULTS_PAGE)
    app.session_state["last_result"] = result
    app.run(timeout=30)
    assert not app.exception

    build_button = next(button for button in app.button if button.label == "Build complete ZIP")
    build_button.click().run(timeout=30)
    assert not app.exception
    assert (output / "complete_run.zip").exists()
