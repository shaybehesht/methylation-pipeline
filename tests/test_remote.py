from pathlib import Path

from app import remote


def test_build_sshfs_command_with_proxy_jump():
    command = remote.build_sshfs_command(
        remote_user="u244415",
        analysis_host="analysis1.hgsc.bcm.edu",
        remote_path="/data/bams",
        mount_point="/home/me/mnt",
        jump_host="login1.hgsc.bcm.edu",
        read_only=True,
    )
    assert command[0] == "sshfs"
    assert command[1] == "u244415@analysis1.hgsc.bcm.edu:/data/bams"
    assert command[2] == "/home/me/mnt"
    text = remote.command_string(command)
    # jump defaults to the same username when jump_user is omitted
    assert "ProxyJump=u244415@login1.hgsc.bcm.edu" in text
    assert "-o ro" in text
    assert "BatchMode" not in text  # interactive by default so password/2FA can be entered


def test_build_sshfs_command_batch_mode_and_jump_user():
    command = remote.build_sshfs_command(
        remote_user="alice", analysis_host="analysis1", remote_path="/x",
        mount_point="/mnt/x", jump_host="login1", jump_user="bob",
        read_only=False, batch_mode=True,
    )
    text = remote.command_string(command)
    assert "ProxyJump=bob@login1" in text
    assert "BatchMode=yes" in text
    assert "-o ro" not in text  # read_only disabled


def test_command_string_quotes_spaces():
    command = remote.build_sshfs_command(
        remote_user="u", analysis_host="h", remote_path="/a b/c",
        mount_point="/mnt/space dir",
    )
    text = remote.command_string(command)
    assert "'u@h:/a b/c'" in text
    assert "'/mnt/space dir'" in text


def test_is_mounted_false_for_plain_dir(tmp_path: Path):
    assert remote.is_mounted(tmp_path) is False
    assert remote.is_mounted(tmp_path / "missing") is False


def test_unmount_command_platform(monkeypatch):
    monkeypatch.setattr(remote.sys, "platform", "darwin")
    assert remote.unmount_command("/mnt/x").startswith("umount ")
    monkeypatch.setattr(remote.sys, "platform", "linux")
    assert remote.unmount_command("/mnt/x").startswith("fusermount -u ")
