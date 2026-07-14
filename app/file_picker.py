"""Server-side local file browser rooted at a configurable data directory.

The picker never uploads files through the browser and never asks the user to
type a path. All navigation is constrained to ``METHYL_TRIO_DATA_ROOT`` (or the
user's home directory when the variable is unset) so a session cannot read
files outside the intended data area.
"""
from __future__ import annotations

import os
from pathlib import Path

import streamlit as st

BAM_EXTENSIONS = (".bam",)
VCF_EXTENSIONS = (".vcf", ".vcf.gz", ".bcf")


def data_root() -> Path:
    """Return the directory that anchors every browse operation.

    Docker deployments mount host data read-only at ``/data`` and set
    ``METHYL_TRIO_DATA_ROOT`` accordingly; native runs default to the user's
    home directory.
    """

    configured = os.environ.get("METHYL_TRIO_DATA_ROOT")
    root = Path(configured).expanduser() if configured else Path.home()
    try:
        return root.resolve()
    except (OSError, RuntimeError):
        return root


def resolve_within(root: str | Path, candidate: str | Path) -> Path | None:
    """Resolve ``candidate`` and return it only if it stays inside ``root``.

    Symlinks are followed before the containment check so a link cannot be used
    to escape the data root. Returns ``None`` when the path escapes the root or
    cannot be resolved.
    """

    root_path = Path(root)
    try:
        resolved_root = root_path.resolve()
        resolved = Path(candidate).resolve()
    except (OSError, RuntimeError):
        return None
    if resolved == resolved_root or resolved_root in resolved.parents:
        return resolved
    return None


def list_entries(
    directory: str | Path, extensions: tuple[str, ...] = ()
) -> tuple[list[Path], list[Path]]:
    """Return sorted (subdirectories, files) inside ``directory``.

    Hidden entries are skipped and files are filtered by ``extensions`` (case
    insensitive) when provided.
    """

    exts = tuple(ext.lower() for ext in extensions)
    directory_path = Path(directory)
    dirs: list[Path] = []
    files: list[Path] = []
    try:
        entries = sorted(directory_path.iterdir(), key=lambda item: item.name.lower())
    except (OSError, PermissionError):
        return dirs, files
    for entry in entries:
        if entry.name.startswith("."):
            continue
        try:
            if entry.is_dir():
                dirs.append(entry)
            elif entry.is_file() and (not exts or entry.name.lower().endswith(exts)):
                files.append(entry)
        except OSError:
            continue
    return dirs, files


def detect_bam_index(bam_path: str | Path) -> Path | None:
    """Return the index for ``bam_path`` if a ``.bam.bai``/``.bai`` file exists.

    ``.csi`` indexes are also recognised. Returns ``None`` when no index is
    present so the caller can explain that one must be created.
    """

    bam = Path(bam_path)
    candidates = [
        bam.with_name(bam.name + ".bai"),
        bam.with_suffix(".bai"),
        bam.with_name(bam.name + ".csi"),
        bam.with_suffix(".csi"),
    ]
    for candidate in candidates:
        try:
            if candidate.is_file():
                return candidate
        except OSError:
            continue
    return None


def _display_path(root: Path, path: Path) -> str:
    try:
        relative = path.relative_to(root)
    except ValueError:
        return "/"
    text = str(relative)
    return "/" if text == "." else f"/{text}"


def _current_dir(dir_state_key: str, root: Path) -> Path:
    stored = st.session_state.get(dir_state_key, str(root))
    current = resolve_within(root, stored)
    if current is None or not current.is_dir():
        current = root
        st.session_state[dir_state_key] = str(root)
    return current


def file_browser(
    label: str,
    *,
    key: str,
    extensions: tuple[str, ...] = (),
    root: str | Path | None = None,
    help: str | None = None,
) -> str:
    """Render a rooted directory browser and return the selected file path.

    The returned value is an empty string until the user selects a file. State
    is namespaced by ``key`` so several independent pickers can coexist on one
    page.
    """

    browse_root = Path(root).resolve() if root is not None else data_root()
    dir_state_key = f"picker_dir::{key}"
    selection_key = f"picker_selected::{key}"
    nav_widget_key = f"picker_nav::{key}"
    file_widget_key = f"picker_file::{key}"

    if dir_state_key not in st.session_state:
        st.session_state[dir_state_key] = str(browse_root)

    current = _current_dir(dir_state_key, browse_root)

    st.markdown(f"**{label}**")
    if help:
        st.caption(help)
    st.caption(f"Browsing: {_display_path(browse_root, current)}")

    subdirs, files = list_entries(current, extensions)

    def _navigate() -> None:
        choice = st.session_state.get(nav_widget_key, "—")
        if not choice or choice == "—":
            return
        here = _current_dir(dir_state_key, browse_root)
        target = here.parent if choice == ".." else here / choice
        resolved = resolve_within(browse_root, target)
        if resolved is not None and resolved.is_dir():
            st.session_state[dir_state_key] = str(resolved)
        st.session_state[nav_widget_key] = "—"

    def _select() -> None:
        choice = st.session_state.get(file_widget_key, "—")
        if not choice or choice == "—":
            return
        here = _current_dir(dir_state_key, browse_root)
        chosen = resolve_within(browse_root, here / choice)
        if chosen is not None and chosen.is_file():
            st.session_state[selection_key] = str(chosen)

    nav_options = ["—"]
    if current != browse_root:
        nav_options.append("..")
    nav_options.extend(entry.name for entry in subdirs)
    st.selectbox(
        "Open folder", nav_options, key=nav_widget_key, on_change=_navigate,
        help="Navigate within the data root. Selecting a folder opens it.",
    )

    file_options = ["—"] + [entry.name for entry in files]
    if len(file_options) > 1:
        st.selectbox(
            "Select file", file_options, key=file_widget_key, on_change=_select,
        )
    else:
        st.caption("No matching files in this folder.")

    selected = st.session_state.get(selection_key, "")
    if selected:
        st.caption(f"Selected: {_display_path(browse_root, Path(selected))}")
    return selected
