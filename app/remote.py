"""Safe helpers for browsing remote data over an SSHFS (FUSE) mount.

Design principle: the application never sees, stores, or transmits a password.
The user authenticates in their own terminal (password + 2FA/Duo stay entirely
between them and the server); this module only *constructs* the exact ``sshfs``
command, verifies that a mount is live, and exposes an optional key/agent-only
mount attempt that can never prompt for or capture a password.

A double SSH hop (``login1`` -> ``analysis1``) is handled with OpenSSH's
``ProxyJump`` so no manual tunnelling is required.
"""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def sshfs_available() -> bool:
    return shutil.which("sshfs") is not None


def _mount_options(read_only: bool, batch_mode: bool, proxy_jump: str | None) -> list[str]:
    options = [
        "reconnect",
        "follow_symlinks",
        "idmap=user",
        "ServerAliveInterval=15",
        "ServerAliveCountMax=3",
    ]
    if read_only:
        options.append("ro")
    if proxy_jump:
        options.append(f"ProxyJump={proxy_jump}")
    if batch_mode:
        # Fail fast instead of ever prompting for a password: only key/agent
        # authentication can succeed with these options.
        options.extend(["BatchMode=yes", "ConnectTimeout=10"])
    return options


def build_sshfs_command(
    *,
    remote_user: str,
    analysis_host: str,
    remote_path: str,
    mount_point: str | Path,
    jump_user: str | None = None,
    jump_host: str | None = None,
    read_only: bool = True,
    batch_mode: bool = False,
) -> list[str]:
    """Return the argv for mounting ``analysis_host:remote_path`` at ``mount_point``.

    When ``jump_host`` is given, the connection is routed through it with
    ``ProxyJump`` (the ``login1`` -> ``analysis1`` hop).
    """

    target = f"{remote_user}@{analysis_host}:{remote_path}"
    proxy_jump = None
    if jump_host:
        proxy_jump = f"{jump_user or remote_user}@{jump_host}"
    command = ["sshfs", target, str(mount_point)]
    for option in _mount_options(read_only, batch_mode, proxy_jump):
        command.extend(["-o", option])
    return command


def command_string(command: list[str]) -> str:
    """Render a command list as a copy-pasteable, safely quoted shell string."""
    return " ".join(shlex.quote(part) for part in command)


def is_mounted(mount_point: str | Path) -> bool:
    path = Path(mount_point)
    try:
        return path.is_dir() and os.path.ismount(str(path))
    except OSError:
        return False


def unmount_command(mount_point: str | Path) -> str:
    quoted = shlex.quote(str(mount_point))
    if sys.platform == "darwin":
        return f"umount {quoted}"
    return f"fusermount -u {quoted}"


def attempt_mount(command: list[str], mount_point: str | Path, timeout: int = 45) -> tuple[bool, str]:
    """Best-effort mount using only key/agent auth (never prompts for a password).

    ``command`` must have been built with ``batch_mode=True`` so it cannot hang
    on an interactive prompt. Returns ``(mounted, output)``.
    """

    Path(mount_point).mkdir(parents=True, exist_ok=True)
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return False, "Timed out. Use the terminal command so you can complete password/2FA login."
    except FileNotFoundError:
        return False, "sshfs is not installed. Install it (see the notes below) or mount manually."
    output = (completed.stdout + completed.stderr).strip()
    return is_mounted(mount_point), output
