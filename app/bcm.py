"""In-app two-hop SSH to BCM HGSC, mirroring the Unicorn connection model.

Connection path: your laptop -> login1.hgsc.bcm.edu -> analysis1.hgsc.bcm.edu

Security model (safe for LOCAL use — each user runs this on their own machine):
  * Each user logs in with their OWN BCM username/password; every command runs
    as that user, so permissions and audit logs stay correct.
  * Credentials live only in Streamlit's in-memory session_state — never written
    to disk, never logged, and cleared on logout.
  * login1's host key is verified against ~/.ssh/known_hosts (rejects MITM). The
    internal login1 -> analysis1 hop is auto-trusted (never leaves BCM's network).

``paramiko`` is imported lazily so importing this module (and the pages that use
it) does not hard-require the dependency until a connection is actually made.
"""
from __future__ import annotations

import posixpath
import stat as stat_module
from dataclasses import dataclass

GATEWAY_HOST = "login1.hgsc.bcm.edu"
TARGET_HOST = "analysis1.hgsc.bcm.edu"
SSH_PORT = 22


@dataclass(frozen=True)
class RemoteEntry:
    name: str
    is_dir: bool
    size: int


def paramiko_available() -> bool:
    try:
        import paramiko  # noqa: F401
    except Exception:
        return False
    return True


def _connect(username: str, password: str, gateway_host: str, target_host: str, port: int = SSH_PORT):
    """Open a two-hop SSH connection; return (gateway, target) clients.

    The gateway client must stay open while the target is used, otherwise the
    tunnel closes underneath it.
    """
    import paramiko

    gateway = paramiko.SSHClient()
    gateway.load_system_host_keys()
    gateway.set_missing_host_key_policy(paramiko.RejectPolicy())
    gateway.connect(
        gateway_host, port=port, username=username, password=password,
        allow_agent=False, look_for_keys=False,
    )
    tunnel = gateway.get_transport().open_channel(
        "direct-tcpip", (target_host, port), ("127.0.0.1", 0),
    )
    target = paramiko.SSHClient()
    target.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    target.connect(
        target_host, port=port, username=username, password=password,
        sock=tunnel, allow_agent=False, look_for_keys=False,
    )
    return gateway, target


def run_command(
    username: str, password: str, command: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> tuple[str, str]:
    gateway, target = _connect(username, password, gateway_host, target_host)
    try:
        _, stdout, stderr = target.exec_command(command)
        return stdout.read().decode(errors="replace"), stderr.read().decode(errors="replace")
    finally:
        target.close()
        gateway.close()


def whoami(username: str, password: str, *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST) -> str:
    out, _ = run_command(
        username, password, "hostname && whoami",
        gateway_host=gateway_host, target_host=target_host,
    )
    return out.strip()


def sort_entries(entries: list[RemoteEntry]) -> list[RemoteEntry]:
    """Directories first, then files, each alphabetical (case-insensitive)."""
    return sorted(entries, key=lambda entry: (not entry.is_dir, entry.name.lower()))


def list_dir(
    username: str, password: str, remote_path: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> list[RemoteEntry]:
    gateway, target = _connect(username, password, gateway_host, target_host)
    try:
        sftp = target.open_sftp()
        entries = [
            RemoteEntry(attr.filename, stat_module.S_ISDIR(attr.st_mode), int(attr.st_size or 0))
            for attr in sftp.listdir_attr(remote_path)
        ]
        sftp.close()
        return sort_entries(entries)
    finally:
        target.close()
        gateway.close()


def exists(
    username: str, password: str, remote_path: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> bool:
    import shlex

    out, _ = run_command(
        username, password, f"test -e {shlex.quote(remote_path)} && echo yes || echo no",
        gateway_host=gateway_host, target_host=target_host,
    )
    return out.strip().endswith("yes")


def bam_index_path(
    username: str, password: str, bam_path: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> str | None:
    """Return the server-side BAM index path (``.bam.bai``/``.bai``/``.csi``) if present."""
    candidates = [bam_path + ".bai", bam_path[:-4] + ".bai" if bam_path.endswith(".bam") else bam_path + ".bai",
                  bam_path + ".csi"]
    for candidate in dict.fromkeys(candidates):
        if exists(username, password, candidate, gateway_host=gateway_host, target_host=target_host):
            return candidate
    return None


def download_to(
    username: str, password: str, remote_file: str, local_path: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> str:
    """Stream a remote file to ``local_path`` (to disk, not memory)."""
    import os

    gateway, target = _connect(username, password, gateway_host, target_host)
    try:
        os.makedirs(os.path.dirname(local_path) or ".", exist_ok=True)
        sftp = target.open_sftp()
        sftp.get(remote_file, local_path)
        sftp.close()
        return local_path
    finally:
        target.close()
        gateway.close()


def join(remote_dir: str, name: str) -> str:
    return posixpath.join(remote_dir, name)


def parent(remote_dir: str) -> str:
    return posixpath.dirname(remote_dir.rstrip("/")) or "/"
