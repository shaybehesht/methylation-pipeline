from pathlib import Path

import streamlit as st

from app import anvil, branding
from app.file_picker import data_roots, register_data_root
from app.state import initialize

initialize()
branding.style()
st.title("🥭 0. AnVIL data (Google Cloud)")
st.caption(
    "Pull modBAMs straight from your AnVIL / Terra workspace bucket on Google "
    "Cloud Storage. Browse the bucket or paste gs:// paths from a Data Table, "
    "localize the trio you need, and it becomes selectable on the Setup page — "
    "all inside AnVIL."
)

info = anvil.workspace_info()

if not anvil.gsutil_available():
    st.warning(
        "`gsutil` was not found on this machine. This page is meant to run inside "
        "an AnVIL / Terra Cloud Environment, where `gsutil` and your Google "
        "credentials are already configured. See `anvil/README.md` for how to "
        "launch MANGO as an interactive app in your workspace."
    )

with st.expander("Workspace and billing", expanded=True):
    st.caption(
        "GREGoR release buckets are **requester-pays**: every read is billed to "
        "your own Terra billing project. Inside a Cloud Environment these fields "
        "are detected automatically; override them only if needed."
    )
    cols = st.columns(2)
    bucket_default = st.session_state.get("anvil_bucket") or info["bucket"] or ""
    project_default = st.session_state.get("anvil_project") or info["project"] or ""
    workspace_bucket = cols[0].text_input(
        "Workspace bucket", bucket_default,
        placeholder="gs://fc-xxxxxxxx-xxxx-xxxx",
        help="The gs:// bucket backing this workspace (WORKSPACE_BUCKET).",
    )
    billing_project = cols[1].text_input(
        "Billing project (requester-pays)", project_default,
        placeholder="my-terra-billing-project",
        help="Google project charged for reads (GOOGLE_PROJECT). Required for "
        "requester-pays consortium buckets.",
    )
    st.session_state.anvil_bucket = workspace_bucket
    st.session_state.anvil_project = billing_project
    if info["name"]:
        st.caption(f"Workspace: {info['namespace'] or '?'}/{info['name']}")

project = billing_project.strip() or None

download_dir = st.text_input(
    "Local download folder",
    st.session_state.get("anvil_download_dir", str(Path.home() / "methyl-trio-anvil")),
    help="Localized BAMs (and their indexes) are saved here, and this folder "
    "becomes browsable in Setup.",
)
st.session_state.anvil_download_dir = download_dir

tab_browse, tab_paste = st.tabs(["Browse bucket", "Paste gs:// paths"])

with tab_browse:
    st.caption(
        "List a bucket or folder, open subfolders, and download the BAMs you need. "
        "Start from the workspace bucket or paste any gs:// prefix."
    )
    start = st.session_state.get("anvil_prefix") or (workspace_bucket.strip() or "")
    location = st.text_input(
        "gs:// location", start,
        placeholder="gs://fc-xxxx/submissions/…  or  gs://bucket/path/to/bams/",
    )
    controls = st.columns(3)
    if controls[0].button("List", type="primary"):
        st.session_state.anvil_prefix = location.strip()
        if not anvil.is_gs_uri(location):
            st.error("Enter a gs:// location to list.")
        else:
            try:
                st.session_state.anvil_listing = anvil.list_objects(
                    location.strip(), project=project
                )
            except Exception as exc:  # noqa: BLE001 - surface gsutil errors
                st.session_state.anvil_listing = []
                st.error(f"Could not list `{location}`: {exc}")
    if controls[1].button("Up one level"):
        current = (st.session_state.get("anvil_prefix") or location).strip().rstrip("/")
        if anvil.is_gs_uri(current) and "/" in current[len("gs://"):]:
            parent = current.rsplit("/", 1)[0] + "/"
            st.session_state.anvil_prefix = parent
            try:
                st.session_state.anvil_listing = anvil.list_objects(parent, project=project)
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not list `{parent}`: {exc}")
        st.rerun()

    filter_text = st.text_input(
        "Filter", key="anvil_filter", placeholder="type to narrow this folder"
    ).strip().lower()

    listing = st.session_state.get("anvil_listing", [])
    if filter_text:
        listing = [entry for entry in listing if filter_text in entry.name.lower()]
    for entry in listing:
        row = st.columns([6, 2])
        icon = "📁" if entry.is_dir else "📄"
        row[0].write(f"{icon} {entry.name}")
        if entry.is_dir:
            if row[1].button("Open", key=f"anvil_open::{entry.uri}"):
                st.session_state.anvil_prefix = entry.uri
                try:
                    st.session_state.anvil_listing = anvil.list_objects(
                        entry.uri, project=project
                    )
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Could not open `{entry.uri}`: {exc}")
                st.rerun()
        elif entry.name.lower().endswith(anvil.BAM_SUFFIXES):
            if row[1].button("Download", key=f"anvil_dl::{entry.uri}"):
                try:
                    with st.spinner(f"Localizing {entry.name}…"):
                        result = anvil.download(entry.uri, download_dir, project=project)
                    register_data_root(download_dir)
                    if result.index:
                        st.success(
                            f"Downloaded {entry.name} and its index. "
                            f"'{download_dir}' is now browsable in Setup."
                        )
                    else:
                        st.warning(
                            f"Downloaded {entry.name}, but no index (.bai/.csi) was "
                            "found beside it. Create one with `samtools index` so "
                            "modkit can read regions efficiently."
                        )
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Download failed: {exc}")

with tab_paste:
    st.caption(
        "Paste the gs:// path of each modBAM (one per line) — for example the "
        "values from a workspace Data Table column — and fetch them all at once. "
        "Matching .bai/.csi indexes are fetched automatically when present."
    )
    uris_text = st.text_area(
        "gs:// BAM paths", st.session_state.get("anvil_uris_text", ""),
        placeholder="gs://fc-xxxx/proband.bam\ngs://fc-xxxx/mother.bam\ngs://fc-xxxx/father.bam",
        height=120,
    )
    st.session_state.anvil_uris_text = uris_text
    if st.button("Fetch all", type="primary", key="anvil_fetch_all"):
        uris = [line.strip() for line in uris_text.splitlines() if line.strip()]
        invalid = [uri for uri in uris if not anvil.is_gs_uri(uri)]
        if not uris:
            st.error("Paste at least one gs:// BAM path.")
        elif invalid:
            st.error("These are not gs:// paths: " + ", ".join(invalid))
        else:
            ok, failed = 0, []
            with st.spinner(f"Localizing {len(uris)} file(s)…"):
                for uri in uris:
                    try:
                        anvil.download(uri, download_dir, project=project)
                        ok += 1
                    except Exception as exc:  # noqa: BLE001
                        failed.append(f"{anvil.basename(uri)}: {exc}")
            if ok:
                register_data_root(download_dir)
                st.success(
                    f"Localized {ok} file(s) to '{download_dir}'. It is now "
                    "browsable in Setup."
                )
            for message in failed:
                st.error(f"Failed — {message}")

st.divider()
st.caption("Currently browsable locations: " + ", ".join(str(root) for root in data_roots()))
st.info(
    "Targeted gene-panel runs only need small region slices, so localizing whole "
    "BAMs is best for a specific trio. For genome-wide analysis of very large "
    "BAMs, consider the WDL workflow path (see `wdl/`) which runs on Cromwell "
    "without downloading BAMs into the app. Write run outputs to a folder inside "
    "the workspace bucket (mounted under your Cloud Environment) so results "
    "persist in AnVIL."
)
