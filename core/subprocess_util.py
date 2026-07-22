"""Shared subprocess helpers for external CLI tools."""
from __future__ import annotations

import shlex
import subprocess
from pathlib import Path


def run_checked(
    command: list[str],
    *,
    log=None,
    log_filepath: str | Path | None = None,
    tool: str,
) -> subprocess.CompletedProcess[str]:
    """Run *command*, append stdout/stderr (and optional log file) to *log*."""
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if log is not None:
        if completed.stdout:
            log.write(completed.stdout)
        if completed.stderr:
            log.write(completed.stderr)
    if completed.returncode:
        parts = [
            f"{tool} failed (exit {completed.returncode})",
            f"command: {shlex.join(command)}",
        ]
        detail = (completed.stderr or completed.stdout or "").strip()
        if detail:
            parts.append("output:\n" + "\n".join(detail.splitlines()[-30:]))
        if log_filepath:
            log_path = Path(log_filepath)
            if log_path.exists():
                log_text = log_path.read_text(encoding="utf-8", errors="replace").strip()
                if log_text:
                    parts.append("log file:\n" + "\n".join(log_text.splitlines()[-30:]))
        raise RuntimeError("\n".join(parts))
    return completed
