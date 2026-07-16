from app import bcm


def test_default_hosts_match_bcm():
    assert bcm.GATEWAY_HOST == "login1.hgsc.bcm.edu"
    assert bcm.TARGET_HOST == "analysis1.hgsc.bcm.edu"
    assert bcm.SSH_PORT == 22


def test_sort_entries_dirs_first_then_alpha():
    entries = [
        bcm.RemoteEntry("zeta.bam", False, 10),
        bcm.RemoteEntry("Alpha", True, 0),
        bcm.RemoteEntry("beta.bam", False, 5),
        bcm.RemoteEntry("gamma", True, 0),
    ]
    ordered = [entry.name for entry in bcm.sort_entries(entries)]
    assert ordered == ["Alpha", "gamma", "beta.bam", "zeta.bam"]


def test_join_and_parent_are_posix():
    assert bcm.join("/stornext/next-gen", "BH16732") == "/stornext/next-gen/BH16732"
    assert bcm.parent("/stornext/next-gen/BH16732") == "/stornext/next-gen"
    assert bcm.parent("/") == "/"


def test_paramiko_available():
    # paramiko is a declared dependency, so it should import in the test env.
    assert bcm.paramiko_available() is True


def test_build_view_command_is_read_only_and_region_scoped():
    cmd = bcm.build_view_command("/stornext/x/proband.bam", ["chr3:100-200", "chrX:5-9"])
    # runs in a login shell so module/conda init is sourced
    assert cmd.startswith("bash -lc ")
    assert "samtools view -b" in cmd
    assert "/stornext/x/proband.bam" in cmd
    assert "chr3:100-200" in cmd and "chrX:5-9" in cmd
    # no write/redirect on the server
    assert ">" not in cmd and " -o " not in cmd


def test_build_view_command_without_login_shell_is_plain():
    cmd = bcm.build_view_command("/x.bam", ["chr1:1-2"], login_shell=False)
    assert cmd.startswith("samtools view -b ")


def test_build_view_command_with_setup_module_load():
    cmd = bcm.build_view_command(
        "/x.bam", ["chr1:1-2"], setup="module load samtools", login_shell=False
    )
    assert cmd.startswith("module load samtools && samtools view -b ")


def test_build_view_command_requires_regions():
    import pytest

    with pytest.raises(ValueError):
        bcm.build_view_command("/x.bam", [])


def test_build_view_command_quotes_spaces():
    cmd = bcm.build_view_command(
        "/a b/proband.bam", ["chr1:1-2"], samtools="/opt/s t/samtools", login_shell=False
    )
    assert "'/a b/proband.bam'" in cmd
    assert "'/opt/s t/samtools'" in cmd


def test_build_pysam_slice_command():
    cmd = bcm.build_pysam_slice_command(
        "python3", "/home/u/.methyl_trio/remote_slice.py",
        "/stornext/x/proband.bam", ["chr3:101-200", "chrX:1-50"],
    )
    assert cmd.startswith("python3 ")
    assert "/home/u/.methyl_trio/remote_slice.py" in cmd
    assert "/stornext/x/proband.bam" in cmd
    assert "chr3:101-200" in cmd and "chrX:1-50" in cmd
    assert ">" not in cmd  # streams to stdout, no server-side write


def test_build_pysam_slice_command_requires_regions():
    import pytest

    with pytest.raises(ValueError):
        bcm.build_pysam_slice_command("python3", "/s.py", "/x.bam", [])


def test_slice_script_is_valid_python():
    import ast

    ast.parse(bcm.SLICE_SCRIPT)
    assert "pysam" in bcm.SLICE_SCRIPT
