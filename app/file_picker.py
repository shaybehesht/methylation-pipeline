"""Reusable server-side local file browser rooted at the data directory."""
from __future__ import annotations

import os
from pathlib import Path

import streamlit as st


def data_root() -> Path:
    """Return the directory users may browse.

    In Docker the host data drive is mounted at ``/data`` and exported through
    ``METHYL_TRIO_DATA_ROOT``. Native runs fall back to the user's home so the
    picker never exposes the whole filesystem.
    """
    configured = os.environ.get("METHYL_TRIO_DATA_ROOT")
    root = Path(configured) if configured else Path.home()
    return root.expanduser().resolve()


def _within_root(candidate: Path, root: Path) -> bool:
    try:
        candidate.resolve().relative_to(root)
        return True
    except ValueError:
        return False


def resolve_in_root(relative: str) -> Path | None:
    """Resolve a user-visible relative path back to an absolute path in root."""
    if not relative:
        return None
    root = data_root()
    candidate = (root / relative).resolve()
    if not _within_root(candidate, root) or not candidate.exists():
        return None
    return candidate


def _display(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root)) or "."
    except ValueError:
        return str(path)


def pick_file(label: str, key: str, extensions: tuple[str, ...], help: str | None = None) -> str | None:
    """Render a rooted browser and return the chosen absolute path as a string.

    Selection state is stored per ``key`` so multiple pickers coexist. Only
    files whose names end with one of ``extensions`` are selectable; directory
    navigation is confined to the data root.
    """
    root = data_root()
    st.markdown(f"**{label}**")
    if help:
        st.caption(help)
    if not root.exists():
        st.error(f"Data directory does not exist: {root}")
        return st.session_state.get(f"picker_selected_{key}")

    cwd_key = f"picker_cwd_{key}"
    selected_key = f"picker_selected_{key}"
    current = Path(st.session_state.get(cwd_key, str(root))).resolve()
    if not _within_root(current, root) and current != root:
        current = root

    st.caption(f"Location: {_display(current, root)}")
    columns = st.columns([1, 4])
    if columns[0].button("Up", key=f"picker_up_{key}", disabled=current == root):
        st.session_state[cwd_key] = str(current.parent if current != root else root)
        st.rerun()

    try:
        entries = sorted(
            current.iterdir(), key=lambda item: (item.is_file(), item.name.lower())
        )
    except PermissionError:
        st.error("Permission denied reading this directory.")
        entries = []

    directories = [item for item in entries if item.is_dir()]
    files = [item for item in entries if item.is_file() and item.name.lower().endswith(extensions)]

    if directories:
        folder_names = ["-"] + [item.name for item in directories]
        chosen = columns[1].selectbox(
            "Open folder", folder_names, key=f"picker_dir_{key}"
        )
        if chosen != "-":
            st.session_state[cwd_key] = str((current / chosen).resolve())
            st.rerun()

    if not files:
        st.info(f"No {', '.join(extensions)} files in this folder.")
    else:
        file_names = ["-"] + [item.name for item in files]
        chosen_file = st.selectbox("Select file", file_names, key=f"picker_file_{key}")
        if chosen_file != "-":
            st.session_state[selected_key] = str((current / chosen_file).resolve())

    selected = st.session_state.get(selected_key)
    if selected and _within_root(Path(selected), root) and Path(selected).exists():
        st.success(f"Selected: {_display(Path(selected), root)}")
        return selected
    return None


def bam_index(bam_path: str) -> str | None:
    """Return the companion index for a BAM if present (.bam.bai or .bai)."""
    bam = Path(bam_path)
    for candidate in (Path(f"{bam}.bai"), bam.with_suffix(".bai")):
        if candidate.exists():
            return str(candidate)
    return None
