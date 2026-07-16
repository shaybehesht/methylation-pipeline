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

# Parents under which operating systems expose removable/external drives. When
# no data root is configured these are auto-detected so external hard drives are
# browsable without any extra setup.
_EXTERNAL_MOUNT_PARENTS = ("/Volumes", "/media", "/run/media", "/mnt")


def _resolve(path: str | Path) -> Path:
    candidate = Path(path).expanduser()
    try:
        return candidate.resolve()
    except (OSError, RuntimeError):
        return candidate


def _external_mount_roots() -> list[Path]:
    """Return existing, non-empty external-mount parents (e.g. ``/Volumes``)."""

    roots: list[Path] = []
    for parent in _EXTERNAL_MOUNT_PARENTS:
        parent_path = Path(parent)
        try:
            if not parent_path.is_dir():
                continue
            has_entries = any(True for _ in parent_path.iterdir())
        except (OSError, PermissionError):
            continue
        if not has_entries:
            continue
        resolved = _resolve(parent_path)
        if resolved not in roots:
            roots.append(resolved)
    return roots


def data_roots() -> list[Path]:
    """Return every directory the file browser is allowed to start from.

    ``METHYL_TRIO_DATA_ROOT`` may list several locations separated by the OS
    path separator (``:`` on POSIX, ``;`` on Windows) — Docker deployments set
    it to ``/data``. When it is unset, browsing starts at the user's home
    directory and any mounted external drives (``/Volumes``, ``/media``,
    ``/run/media``, ``/mnt``) are added automatically.
    """

    configured = os.environ.get("METHYL_TRIO_DATA_ROOT")
    roots: list[Path] = []
    if configured:
        for part in configured.split(os.pathsep):
            part = part.strip()
            if not part:
                continue
            resolved = _resolve(part)
            if resolved not in roots:
                roots.append(resolved)
    else:
        roots.append(_resolve(Path.home()))
        for mount in _external_mount_roots():
            if mount not in roots:
                roots.append(mount)
    for extra in _session_roots():
        if extra not in roots and extra.exists():
            roots.append(extra)
    return roots or [_resolve(Path.home())]


_EXTRA_ROOTS_KEY = "methyl_trio_extra_roots"


def _session_roots() -> list[Path]:
    """Return roots registered at runtime (e.g. mounted remote data)."""
    try:
        registered = st.session_state.get(_EXTRA_ROOTS_KEY, [])
    except Exception:
        return []
    return [_resolve(path) for path in registered]


def register_data_root(path: str | Path) -> Path:
    """Register an extra browsable root for this session (e.g. an SSHFS mount)."""
    resolved = _resolve(path)
    try:
        registered = st.session_state.setdefault(_EXTRA_ROOTS_KEY, [])
        if str(resolved) not in registered:
            registered.append(str(resolved))
    except Exception:
        pass
    return resolved


def data_root() -> Path:
    """Return the primary browse root (the first configured data root)."""

    return data_roots()[0]


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


def _normalize_roots(roots: object) -> list[Path]:
    if roots is None:
        return data_roots()
    if isinstance(roots, (str, Path)):
        roots = [roots]
    resolved: list[Path] = []
    for item in roots:
        candidate = _resolve(item)
        if candidate not in resolved:
            resolved.append(candidate)
    return resolved or data_roots()


def file_browser(
    label: str,
    *,
    key: str,
    extensions: tuple[str, ...] = (),
    roots: object = None,
    help: str | None = None,
) -> str:
    """Render a rooted directory browser and return the selected file path.

    Browsing is confined to one of the configured data roots (home directory,
    ``METHYL_TRIO_DATA_ROOT`` entries, and detected external drives). When more
    than one root is available a "Location" selector chooses between them. The
    returned value is an empty string until the user selects a file. State is
    namespaced by ``key`` so several independent pickers can coexist on one
    page.
    """

    root_list = _normalize_roots(roots)
    dir_state_key = f"picker_dir::{key}"
    selection_key = f"picker_selected::{key}"
    nav_widget_key = f"picker_nav::{key}"
    file_widget_key = f"picker_file::{key}"
    root_widget_key = f"picker_root::{key}"

    st.markdown(f"**{label}**")
    if help:
        st.caption(help)

    if len(root_list) > 1:
        root_options = [str(root) for root in root_list]
        if st.session_state.get(root_widget_key) not in root_options:
            st.session_state[root_widget_key] = root_options[0]

        def _switch_root() -> None:
            st.session_state[dir_state_key] = st.session_state[root_widget_key]
            st.session_state[nav_widget_key] = "—"
            st.session_state.pop(file_widget_key, None)

        st.selectbox(
            "Location", root_options, key=root_widget_key, on_change=_switch_root,
            help="Choose which data root or external drive to browse.",
        )
        browse_root = _resolve(st.session_state[root_widget_key])
    else:
        browse_root = root_list[0]

    if dir_state_key not in st.session_state:
        st.session_state[dir_state_key] = str(browse_root)

    current = _current_dir(dir_state_key, browse_root)

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
