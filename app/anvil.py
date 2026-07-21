"""Read modBAMs from an AnVIL / Terra workspace on Google Cloud Storage.

MANGO runs unchanged inside an AnVIL interactive Cloud Environment. The data in
a GREGoR (or any AnVIL) workspace lives in a Google Cloud Storage bucket and is
referenced by ``gs://`` URIs in the workspace Data Tables. This module lets the
app browse that bucket and localize the specific BAMs (and their indexes) a trio
needs into a folder that then becomes browsable on the Setup page — mirroring the
"Remote data" (SSH) flow but for Google Cloud.

Design notes:

* All Google access goes through the ``gsutil``/``gcloud`` command-line tools,
  which are pre-installed and pre-authenticated inside every Terra Cloud
  Environment. No extra Python client library or service-account key is needed,
  and nothing is uploaded through the browser.
* GREGoR consortium release buckets are **requester-pays**, so every read is
  billed to the user's own Terra billing project. That project id is passed with
  ``gsutil -u <project>`` (it is available as ``GOOGLE_PROJECT`` in the
  environment). Without it, listing or copying a requester-pays object fails.
* The command builders and parsers are pure functions so they can be unit
  tested without a network or the gcloud SDK; the functions that actually shell
  out accept an injectable ``run`` callable for the same reason.
"""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

# A runner takes an argv list and returns (returncode, stdout, stderr).
Runner = Callable[[list[str]], "tuple[int, str, str]"]

# Terra/AnVIL Cloud Environments export these describing the current workspace.
ENV_BUCKET = "WORKSPACE_BUCKET"
ENV_PROJECT = "GOOGLE_PROJECT"
ENV_NAMESPACE = "WORKSPACE_NAMESPACE"
ENV_NAME = "WORKSPACE_NAME"
ENV_OWNER = "OWNER_EMAIL"

BAM_SUFFIXES = (".bam",)
# Index siblings tried, in order of preference, next to a ``*.bam`` object.
_INDEX_SUFFIXES = (".bam.bai", ".bai", ".bam.csi", ".csi")


@dataclass(frozen=True)
class GcsEntry:
    """A single object or prefix returned by listing a ``gs://`` location."""

    name: str
    uri: str
    is_dir: bool


@dataclass(frozen=True)
class DownloadResult:
    """Local paths produced by localizing one BAM (and any index found)."""

    bam: str
    index: str | None
    index_uri: str | None


def _run(command: list[str], timeout: int | None = None) -> tuple[int, str, str]:
    """Default runner: execute ``command`` and capture its output."""
    completed = subprocess.run(
        command, capture_output=True, text=True, timeout=timeout, check=False
    )
    return completed.returncode, completed.stdout, completed.stderr


def gsutil_available() -> bool:
    """True when the ``gsutil`` CLI is on PATH (always so inside Terra)."""
    from shutil import which

    return which("gsutil") is not None


def is_gs_uri(text: str) -> bool:
    return isinstance(text, str) and text.strip().startswith("gs://")


def parse_gs_uri(uri: str) -> tuple[str, str]:
    """Split ``gs://bucket/key`` into ``(bucket, key)``.

    ``key`` may be empty (bucket root) and keeps any trailing slash. Raises
    ``ValueError`` for anything that is not a ``gs://`` URI.
    """
    text = (uri or "").strip()
    if not text.startswith("gs://"):
        raise ValueError(f"Not a gs:// URI: {uri!r}")
    remainder = text[len("gs://"):]
    if not remainder or remainder.startswith("/"):
        raise ValueError(f"gs:// URI is missing a bucket name: {uri!r}")
    bucket, _, key = remainder.partition("/")
    return bucket, key


def basename(uri: str) -> str:
    """Return the final path component of a ``gs://`` URI (no trailing slash)."""
    return uri.rstrip("/").rsplit("/", 1)[-1]


def workspace_bucket() -> str | None:
    """The ``gs://fc-...`` bucket of the current Terra workspace, if known."""
    value = os.environ.get(ENV_BUCKET, "").strip()
    return value or None


def billing_project() -> str | None:
    """The Google project used for requester-pays billing, if known."""
    for key in (ENV_PROJECT, ENV_NAMESPACE):
        value = os.environ.get(key, "").strip()
        if value:
            return value
    return None


def workspace_info() -> dict[str, str | None]:
    """Return the workspace identity exported by the Terra environment."""
    return {
        "bucket": workspace_bucket(),
        "project": billing_project(),
        "namespace": os.environ.get(ENV_NAMESPACE, "").strip() or None,
        "name": os.environ.get(ENV_NAME, "").strip() or None,
        "owner": os.environ.get(ENV_OWNER, "").strip() or None,
    }


def _user_project_flag(project: str | None) -> list[str]:
    """The ``-u <project>`` flag for requester-pays access, or nothing."""
    project = (project or "").strip()
    return ["-u", project] if project else []


def build_ls_command(uri: str, project: str | None = None) -> list[str]:
    """Build the ``gsutil ls`` argv for ``uri`` (with requester-pays billing)."""
    return ["gsutil", *_user_project_flag(project), "ls", uri]


def build_cp_command(source: str, destination: str, project: str | None = None) -> list[str]:
    """Build the ``gsutil cp`` argv (with requester-pays billing)."""
    return ["gsutil", *_user_project_flag(project), "cp", source, destination]


def build_stat_command(uri: str, project: str | None = None) -> list[str]:
    """Build a quiet ``gsutil stat`` argv used to test whether an object exists."""
    return ["gsutil", *_user_project_flag(project), "-q", "stat", uri]


def index_uri_candidates(bam_uri: str) -> list[str]:
    """Return the candidate index URIs to look for next to ``bam_uri``.

    ``sample.bam`` yields ``sample.bam.bai``, ``sample.bai``, ``sample.bam.csi``
    and ``sample.csi`` in preference order.
    """
    text = bam_uri.strip()
    lower = text.lower()
    stem = text[:-4] if lower.endswith(".bam") else text
    return [
        f"{text}.bai",
        f"{stem}.bai",
        f"{text}.csi",
        f"{stem}.csi",
    ]


def parse_ls_output(stdout: str, listed_uri: str = "") -> list[GcsEntry]:
    """Turn ``gsutil ls`` stdout into sorted directory + file entries.

    Lines ending in ``/`` are prefixes (folders); everything else is an object.
    The ``listed_uri`` itself is skipped so a self-referential line does not
    appear as a child entry.
    """
    listed = listed_uri.strip().rstrip("/")
    dirs: list[GcsEntry] = []
    files: list[GcsEntry] = []
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line or not line.startswith("gs://"):
            continue
        if line.rstrip("/") == listed:
            continue
        is_dir = line.endswith("/")
        name = basename(line) + ("/" if is_dir else "")
        entry = GcsEntry(name=name, uri=line, is_dir=is_dir)
        (dirs if is_dir else files).append(entry)
    dirs.sort(key=lambda item: item.name.lower())
    files.sort(key=lambda item: item.name.lower())
    return dirs + files


def list_objects(
    uri: str, project: str | None = None, run: Runner = _run
) -> list[GcsEntry]:
    """List the folders and objects directly under ``uri``.

    Raises ``RuntimeError`` with gsutil's message on failure (for example a
    requester-pays bucket accessed without a billing project).
    """
    code, out, err = run(build_ls_command(uri, project))
    if code != 0:
        raise RuntimeError((err or out or f"gsutil ls failed for {uri}").strip())
    return parse_ls_output(out, listed_uri=uri)


def object_exists(uri: str, project: str | None = None, run: Runner = _run) -> bool:
    """True when ``uri`` is an existing object (via ``gsutil stat``)."""
    code, _out, _err = run(build_stat_command(uri, project))
    return code == 0


def first_existing_index(
    bam_uri: str, project: str | None = None, run: Runner = _run
) -> str | None:
    """Return the first existing index sibling of ``bam_uri`` (or ``None``)."""
    for candidate in index_uri_candidates(bam_uri):
        if object_exists(candidate, project, run):
            return candidate
    return None


def download(
    bam_uri: str,
    destination_dir: str | Path,
    project: str | None = None,
    with_index: bool = True,
    run: Runner = _run,
) -> DownloadResult:
    """Copy ``bam_uri`` (and its index, if present) into ``destination_dir``.

    The BAM copy failing raises ``RuntimeError``; a missing or un-copyable index
    is tolerated (the Setup page then explains how to build one) so a genome
    without a co-located index still localizes.
    """
    dest_dir = Path(destination_dir).expanduser()
    dest_dir.mkdir(parents=True, exist_ok=True)
    local_bam = dest_dir / basename(bam_uri)

    code, out, err = run(build_cp_command(bam_uri, str(local_bam), project))
    if code != 0:
        raise RuntimeError((err or out or f"gsutil cp failed for {bam_uri}").strip())

    index_local: str | None = None
    index_uri: str | None = None
    if with_index:
        index_uri = first_existing_index(bam_uri, project, run)
        if index_uri:
            local_index = dest_dir / basename(index_uri)
            code, _out, _err = run(build_cp_command(index_uri, str(local_index), project))
            if code == 0:
                index_local = str(local_index)
            else:
                index_uri = None
    return DownloadResult(bam=str(local_bam), index=index_local, index_uri=index_uri)
