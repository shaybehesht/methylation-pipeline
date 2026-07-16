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


def build_view_command(
    remote_bam: str, regions: list[str], samtools: str = "samtools",
    *, setup: str | None = None, login_shell: bool = True,
) -> str:
    """Build a read-only ``samtools view -b`` command for specific regions.

    Regions are required so the whole BAM is never streamed by accident. The
    command only reads on the server; its BAM output goes to stdout (streamed to
    a local file), so no write permission on the server is needed.

    ``setup`` runs first (e.g. ``module load samtools`` or ``source
    ~/.bashrc``). ``login_shell`` wraps everything in ``bash -lc`` so the
    server's profile/module initialization is sourced — this fixes
    ``samtools: command not found`` under non-interactive SSH.
    """
    import shlex

    if not regions:
        raise ValueError("At least one region is required to slice a BAM")
    region_args = " ".join(shlex.quote(region) for region in regions)
    inner = f"{shlex.quote(samtools)} view -b {shlex.quote(remote_bam)} {region_args}"
    if setup and setup.strip():
        inner = f"{setup.strip()} && {inner}"
    if login_shell:
        return f"bash -lc {shlex.quote(inner)}"
    return inner


def slice_bam(
    username: str, password: str, remote_bam: str, regions: list[str], local_bam: str,
    *, samtools: str = "samtools", setup: str | None = None, login_shell: bool = True,
    gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST, chunk: int = 1 << 16,
) -> str:
    """Stream a region-restricted BAM from the server to ``local_bam`` and index it.

    Reads only on the server (no writes there); transfers only the requested
    regions. The local file is coordinate-sorted (it inherits the server BAM's
    order) and indexed with pysam.
    """
    import os

    import pysam

    command = build_view_command(
        remote_bam, regions, samtools=samtools, setup=setup, login_shell=login_shell
    )
    gateway, target = _connect(username, password, gateway_host, target_host)
    try:
        _, stdout, stderr = target.exec_command(command)
        os.makedirs(os.path.dirname(local_bam) or ".", exist_ok=True)
        with open(local_bam, "wb") as handle:
            while True:
                data = stdout.read(chunk)
                if not data:
                    break
                handle.write(data)
        error = stderr.read().decode(errors="replace")
        status = stdout.channel.recv_exit_status()
        if status != 0:
            raise RuntimeError(error.strip() or f"samtools exited with code {status}")
    finally:
        target.close()
        gateway.close()
    pysam.index(str(local_bam))
    return local_bam


def locate_tool(
    username: str, password: str, tool: str = "samtools",
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> str:
    """Return discovery output to help find a tool on the server."""
    import shlex

    quoted = shlex.quote(tool)
    discovery = (
        f"command -v {quoted} || true; "
        f"module avail 2>&1 | grep -i {quoted} || true; "
        f"ls /opt/*/bin/{tool} /usr/local/bin/{tool} 2>/dev/null || true"
    )
    out, err = run_command(
        username, password, f"bash -lc {shlex.quote(discovery)}",
        gateway_host=gateway_host, target_host=target_host,
    )
    return (out + ("\n" + err if err.strip() else "")).strip() or "Nothing found."


REMOTE_APP_DIR = ".methyl_trio"
REMOTE_SLICER_NAME = "remote_slice.py"

# Runs on the server with only python3 + pysam (installed into the user's home).
# Streams a region-restricted, coordinate-sorted BAM to stdout — no server-side
# write, no samtools required.
SLICE_SCRIPT = '''#!/usr/bin/env python3
"""Stream a region-restricted BAM to stdout using pysam (no samtools needed)."""
import sys

import pysam


def main() -> int:
    if len(sys.argv) < 3:
        sys.stderr.write("usage: remote_slice.py BAM REGION [REGION ...]\\n")
        return 2
    bam_path = sys.argv[1]
    regions = sys.argv[2:]
    infile = pysam.AlignmentFile(bam_path, "rb")
    out = pysam.AlignmentFile("-", "wb", template=infile)
    seen = set()
    for region in regions:
        for read in infile.fetch(region=region):
            key = (read.query_name, read.reference_id, read.reference_start, read.flag)
            if key in seen:
                continue
            seen.add(key)
            out.write(read)
    out.close()
    infile.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''


def home_dir(
    username: str, password: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> str:
    out, _ = run_command(
        username, password, "echo $HOME",
        gateway_host=gateway_host, target_host=target_host,
    )
    return out.strip() or f"/home/{username}"


def _run_on_gateway(username: str, password: str, command: str, *, gateway_host: str, port: int = SSH_PORT):
    import paramiko

    gateway = paramiko.SSHClient()
    gateway.load_system_host_keys()
    gateway.set_missing_host_key_policy(paramiko.RejectPolicy())
    gateway.connect(
        gateway_host, port=port, username=username, password=password,
        allow_agent=False, look_for_keys=False,
    )
    try:
        _, stdout, stderr = gateway.exec_command(command)
        return stdout.read().decode(errors="replace"), stderr.read().decode(errors="replace")
    finally:
        gateway.close()


def ensure_remote_slicer(
    username: str, password: str, *, python: str = "python3",
    gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> tuple[bool, str]:
    """Upload the pysam slicer to the user's home and make sure pysam imports.

    pysam is installed with ``pip install --user`` — tried on the analysis host
    first, then on the gateway (login1), which shares the home filesystem and
    usually has outbound internet. Returns ``(ready, log)``.
    """
    import shlex

    log: list[str] = []
    home = home_dir(username, password, gateway_host=gateway_host, target_host=target_host)
    remote_dir = f"{home}/{REMOTE_APP_DIR}"
    script_path = f"{remote_dir}/{REMOTE_SLICER_NAME}"

    gateway, target = _connect(username, password, gateway_host, target_host)
    try:
        target.exec_command(f"mkdir -p {shlex.quote(remote_dir)}")[1].channel.recv_exit_status()
        sftp = target.open_sftp()
        with sftp.open(script_path, "w") as handle:
            handle.write(SLICE_SCRIPT)
        sftp.close()
    finally:
        target.close()
        gateway.close()
    log.append(f"Uploaded slicer to {script_path}")

    check = f"bash -lc {shlex.quote(python + ' -c \"import pysam\" && echo PYSAM_OK')}"

    def has_pysam() -> bool:
        out, _ = run_command(username, password, check, gateway_host=gateway_host, target_host=target_host)
        return "PYSAM_OK" in out

    if has_pysam():
        log.append("pysam already available on the analysis host.")
        return True, "\n".join(log)

    install = f"bash -lc {shlex.quote(python + ' -m pip install --user --quiet pysam')}"
    log.append("Installing pysam --user on the analysis host…")
    run_command(username, password, install, gateway_host=gateway_host, target_host=target_host)
    if has_pysam():
        log.append("pysam installed on the analysis host.")
        return True, "\n".join(log)

    log.append("Analysis host could not install pysam (likely no internet); trying the login host…")
    try:
        _run_on_gateway(username, password, install, gateway_host=gateway_host)
    except Exception as exc:  # noqa: BLE001
        log.append(f"Gateway install error: {exc}")
    if has_pysam():
        log.append("pysam installed via the login host (shared home).")
        return True, "\n".join(log)

    log.append(
        "Could not make pysam importable automatically. On a host with internet "
        "and your shared home, run:  python3 -m pip install --user pysam"
    )
    return False, "\n".join(log)


def build_pysam_slice_command(python: str, script_path: str, remote_bam: str, regions: list[str]) -> str:
    import shlex

    if not regions:
        raise ValueError("At least one region is required to slice a BAM")
    region_args = " ".join(shlex.quote(region) for region in regions)
    return (
        f"{shlex.quote(python)} {shlex.quote(script_path)} "
        f"{shlex.quote(remote_bam)} {region_args}"
    )


def slice_bam_pysam(
    username: str, password: str, remote_bam: str, regions: list[str], local_bam: str,
    *, python: str = "python3", gateway_host: str = GATEWAY_HOST,
    target_host: str = TARGET_HOST, chunk: int = 1 << 16,
) -> str:
    """Stream a region-restricted BAM produced by the remote pysam slicer to disk."""
    import os

    import pysam

    home = home_dir(username, password, gateway_host=gateway_host, target_host=target_host)
    script_path = f"{home}/{REMOTE_APP_DIR}/{REMOTE_SLICER_NAME}"
    command = build_pysam_slice_command(python, script_path, remote_bam, regions)
    gateway, target = _connect(username, password, gateway_host, target_host)
    try:
        _, stdout, stderr = target.exec_command(command)
        os.makedirs(os.path.dirname(local_bam) or ".", exist_ok=True)
        with open(local_bam, "wb") as handle:
            while True:
                data = stdout.read(chunk)
                if not data:
                    break
                handle.write(data)
        error = stderr.read().decode(errors="replace")
        status = stdout.channel.recv_exit_status()
        if status != 0:
            raise RuntimeError(error.strip() or f"remote slicer exited with code {status}")
    finally:
        target.close()
        gateway.close()
    pysam.index(str(local_bam))
    return local_bam


def diagnose_server(
    username: str, password: str,
    *, gateway_host: str = GATEWAY_HOST, target_host: str = TARGET_HOST,
) -> str:
    """Report which analysis tools exist, conda envs, and whether HOME is writable."""
    import shlex

    script = r"""
echo "== tools on PATH ==";
for t in samtools bcftools modkit tabix bgzip python3; do
  printf "%s: " "$t"; command -v "$t" 2>/dev/null || echo "(missing)";
done
echo; echo "== module spider samtools ==";
{ module spider samtools; } 2>&1 | head -20 || true
echo; echo "== conda envs ==";
ls -d ~/miniconda3/envs/* ~/anaconda3/envs/* ~/.conda/envs/* 2>/dev/null || echo "(none)"
echo; echo "== samtools under home/conda ==";
ls ~/*/bin/samtools ~/.conda/envs/*/bin/samtools 2>/dev/null || echo "(none)"
echo; echo "== HOME writable? ==";
{ f="$HOME/.methyl_trio_write_test"; touch "$f" 2>/dev/null && echo "HOME WRITABLE" && rm -f "$f"; } || echo "HOME NOT WRITABLE"
echo; echo "== software trees ==";
ls -d /hgsc_software/* /stornext/*/software 2>/dev/null | head -20 || true
"""
    out, err = run_command(
        username, password, f"bash -lc {shlex.quote(script)}",
        gateway_host=gateway_host, target_host=target_host,
    )
    return (out + ("\n" + err if err.strip() else "")).strip() or "No output."


def join(remote_dir: str, name: str) -> str:
    return posixpath.join(remote_dir, name)


def parent(remote_dir: str) -> str:
    return posixpath.dirname(remote_dir.rstrip("/")) or "/"
