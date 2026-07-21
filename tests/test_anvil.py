import pytest

from app import anvil


def test_is_gs_uri_and_parse():
    assert anvil.is_gs_uri("gs://bucket/x.bam")
    assert not anvil.is_gs_uri("/local/x.bam")
    assert not anvil.is_gs_uri("https://x/y")
    assert anvil.parse_gs_uri("gs://fc-123/path/to/x.bam") == ("fc-123", "path/to/x.bam")
    assert anvil.parse_gs_uri("gs://fc-123") == ("fc-123", "")
    assert anvil.parse_gs_uri("gs://fc-123/sub/") == ("fc-123", "sub/")


@pytest.mark.parametrize("bad", ["/local", "gs://", "gs:///nobucket", "x"])
def test_parse_gs_uri_rejects_bad(bad):
    with pytest.raises(ValueError):
        anvil.parse_gs_uri(bad)


def test_basename():
    assert anvil.basename("gs://b/a/c/proband.bam") == "proband.bam"
    assert anvil.basename("gs://b/folder/") == "folder"


def test_index_uri_candidates():
    assert anvil.index_uri_candidates("gs://b/x/sample.bam") == [
        "gs://b/x/sample.bam.bai",
        "gs://b/x/sample.bai",
        "gs://b/x/sample.bam.csi",
        "gs://b/x/sample.csi",
    ]


def test_command_builders_include_requester_pays_project():
    assert anvil.build_ls_command("gs://b/x/") == ["gsutil", "ls", "gs://b/x/"]
    assert anvil.build_ls_command("gs://b/x/", "proj") == [
        "gsutil", "-u", "proj", "ls", "gs://b/x/",
    ]
    assert anvil.build_cp_command("gs://b/x.bam", "/tmp/x.bam", "proj") == [
        "gsutil", "-u", "proj", "cp", "gs://b/x.bam", "/tmp/x.bam",
    ]
    assert anvil.build_stat_command("gs://b/x.bam", "proj") == [
        "gsutil", "-u", "proj", "-q", "stat", "gs://b/x.bam",
    ]


def test_parse_ls_output_classifies_and_skips_self():
    stdout = (
        "gs://b/data/\n"
        "gs://b/data/sub/\n"
        "gs://b/data/mother.bam\n"
        "gs://b/data/proband.bam\n"
        "\n"
        "not-a-uri\n"
    )
    entries = anvil.parse_ls_output(stdout, listed_uri="gs://b/data/")
    names = [(e.name, e.is_dir) for e in entries]
    # the listed prefix itself is skipped; dirs sort before files
    assert names == [
        ("sub/", True),
        ("mother.bam", False),
        ("proband.bam", False),
    ]


def test_list_objects_uses_runner_and_raises_on_error():
    calls = []

    def ok_runner(cmd):
        calls.append(cmd)
        return 0, "gs://b/x/a.bam\ngs://b/x/sub/\n", ""

    entries = anvil.list_objects("gs://b/x/", project="proj", run=ok_runner)
    assert calls == [["gsutil", "-u", "proj", "ls", "gs://b/x/"]]
    assert {e.name for e in entries} == {"a.bam", "sub/"}

    def bad_runner(cmd):
        return 1, "", "Bucket is requester pays. Provide a billing project."

    with pytest.raises(RuntimeError, match="requester pays"):
        anvil.list_objects("gs://b/x/", run=bad_runner)


def test_first_existing_index_returns_first_hit():
    def runner(cmd):
        # cmd = [gsutil, -q, stat, <uri>]; only the ".bam.bai" sibling exists
        uri = cmd[-1]
        return (0 if uri.endswith("sample.bam.bai") else 1), "", ""

    assert anvil.first_existing_index("gs://b/sample.bam", run=runner) == "gs://b/sample.bam.bai"

    assert anvil.first_existing_index("gs://b/sample.bam", run=lambda c: (1, "", "")) is None


def test_download_fetches_bam_and_index(tmp_path):
    seen = []

    def runner(cmd):
        seen.append(cmd)
        if cmd[1] == "-q" and cmd[2] == "stat":  # index probe
            return (0 if cmd[-1].endswith("proband.bam.bai") else 1), "", ""
        return 0, "", ""  # cp / ls succeed

    result = anvil.download(
        "gs://b/x/proband.bam", tmp_path, project="proj", run=runner
    )
    assert result.bam == str(tmp_path / "proband.bam")
    assert result.index == str(tmp_path / "proband.bam.bai")
    assert result.index_uri == "gs://b/x/proband.bam.bai"
    # the BAM was copied with the billing project
    assert ["gsutil", "-u", "proj", "cp", "gs://b/x/proband.bam", str(tmp_path / "proband.bam")] in seen


def test_download_tolerates_missing_index(tmp_path):
    def runner(cmd):
        if cmd[1] == "-q":  # no index exists
            return 1, "", ""
        return 0, "", ""

    result = anvil.download("gs://b/x/proband.bam", tmp_path, run=runner)
    assert result.bam == str(tmp_path / "proband.bam")
    assert result.index is None
    assert result.index_uri is None


def test_download_raises_when_bam_copy_fails(tmp_path):
    def runner(cmd):
        return 1, "", "AccessDeniedException: 403"

    with pytest.raises(RuntimeError, match="AccessDenied"):
        anvil.download("gs://b/x/proband.bam", tmp_path, run=runner)


def test_workspace_info_reads_environment(monkeypatch):
    monkeypatch.setenv(anvil.ENV_BUCKET, "gs://fc-abc")
    monkeypatch.setenv(anvil.ENV_PROJECT, "billing-proj")
    monkeypatch.setenv(anvil.ENV_NAME, "GREGoR_Release")
    monkeypatch.setenv(anvil.ENV_NAMESPACE, "anvil-namespace")
    assert anvil.workspace_bucket() == "gs://fc-abc"
    assert anvil.billing_project() == "billing-proj"
    info = anvil.workspace_info()
    assert info["bucket"] == "gs://fc-abc"
    assert info["name"] == "GREGoR_Release"
    assert info["namespace"] == "anvil-namespace"


def test_billing_project_falls_back_to_namespace(monkeypatch):
    monkeypatch.delenv(anvil.ENV_PROJECT, raising=False)
    monkeypatch.setenv(anvil.ENV_NAMESPACE, "ns-as-project")
    assert anvil.billing_project() == "ns-as-project"
