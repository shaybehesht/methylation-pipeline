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
