from pathlib import Path

import streamlit as st

from app import bcm, remote
from app.file_picker import data_roots, register_data_root, register_remote_mapping
from app.state import initialize
from core.annotations import extract_to_regions, panel_regions

initialize()
st.title("0. Remote data (BCM / SSH)")
st.caption(
    "Connect to a remote server, browse BAMs by path, and either download the "
    "ones you need or mount the folder for in-place analysis."
)

tab_login, tab_mount = st.tabs(["Connect to BCM (in-app SSH)", "Mount with SSHFS"])

with tab_login:
    st.caption(
        "Log in with your own BCM credentials over a two-hop SSH connection "
        "(login host → analysis host), browse the server, and download the files "
        "you need. Credentials stay in memory only and are cleared on logout — "
        "nothing is written to disk."
    )
    if not bcm.paramiko_available():
        st.warning("`paramiko` is not installed. Run `pip install -e .` in the repo to add it.")

    default_gateway = st.session_state.get("bcm_gateway", bcm.GATEWAY_HOST)
    default_target = st.session_state.get("bcm_target", bcm.TARGET_HOST)

    if not st.session_state.get("bcm_authed"):
        hosts = st.columns(2)
        gateway_host = hosts[0].text_input("Login / jump host", default_gateway)
        target_host = hosts[1].text_input("Analysis / data host", default_target)
        with st.form("bcm_login"):
            username = st.text_input("BCM username")
            password = st.text_input("BCM password", type="password")
            submitted = st.form_submit_button("Connect")
        if submitted:
            try:
                identity = bcm.whoami(
                    username, password, gateway_host=gateway_host, target_host=target_host
                )
                st.session_state.update({
                    "bcm_authed": True, "bcm_user": username, "bcm_pw": password,
                    "bcm_gateway": gateway_host, "bcm_target": target_host,
                    "bcm_cwd": f"/home/{username}",
                })
                st.success(f"Connected: {identity}")
                st.rerun()
            except Exception as exc:  # noqa: BLE001 - surface any auth/connection error
                st.error(f"Login failed: {exc}")
    else:
        user = st.session_state.bcm_user
        pw = st.session_state.bcm_pw
        gateway_host = st.session_state.bcm_gateway
        target_host = st.session_state.bcm_target

        header = st.columns([4, 1])
        header[0].success(f"Logged in as **{user}** on `{target_host}`")
        if header[1].button("Log out"):
            for key in ("bcm_authed", "bcm_user", "bcm_pw", "bcm_cwd", "bcm_listing"):
                st.session_state.pop(key, None)
            st.rerun()

        remote_dir = st.text_input(
            "Remote directory", st.session_state.get("bcm_cwd", f"/home/{user}"),
            help="Paste any path on the analysis host, e.g. "
            "/stornext/snfs190/next-gen/ONT_trio_analysis/…",
        )
        controls = st.columns(3)
        if controls[0].button("List directory"):
            st.session_state.bcm_cwd = remote_dir
            try:
                st.session_state.bcm_listing = bcm.list_dir(
                    user, pw, remote_dir, gateway_host=gateway_host, target_host=target_host
                )
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not list `{remote_dir}`: {exc}")
        if controls[1].button("Up one level"):
            st.session_state.bcm_cwd = bcm.parent(st.session_state.get("bcm_cwd", remote_dir))
            try:
                st.session_state.bcm_listing = bcm.list_dir(
                    user, pw, st.session_state.bcm_cwd,
                    gateway_host=gateway_host, target_host=target_host,
                )
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not list directory: {exc}")

        download_dir = st.text_input(
            "Local download folder",
            st.session_state.get("bcm_download_dir", str(Path.home() / "methyl-trio-downloads")),
            help="Downloaded BAMs (and their indexes) are saved here and this folder "
            "becomes browsable in Setup.",
        )
        st.session_state.bcm_download_dir = download_dir

        cwd = st.session_state.get("bcm_cwd", remote_dir)
        for entry in st.session_state.get("bcm_listing", []):
            row = st.columns([5, 2, 2])
            icon = "📁" if entry.is_dir else "📄"
            row[0].write(f"{icon} {entry.name}")
            row[1].write("" if entry.is_dir else f"{entry.size:,} B")
            full = bcm.join(cwd, entry.name)
            if entry.is_dir:
                if row[2].button("Open", key=f"open::{full}"):
                    st.session_state.bcm_cwd = full
                    try:
                        st.session_state.bcm_listing = bcm.list_dir(
                            user, pw, full, gateway_host=gateway_host, target_host=target_host
                        )
                    except Exception as exc:  # noqa: BLE001
                        st.error(f"Could not open `{full}`: {exc}")
                    st.rerun()
            elif row[2].button("Download", key=f"dl::{full}"):
                target_local = Path(download_dir) / entry.name
                try:
                    with st.spinner(f"Downloading {entry.name}…"):
                        bcm.download_to(
                            user, pw, full, str(target_local),
                            gateway_host=gateway_host, target_host=target_host,
                        )
                        if entry.name.endswith(".bam"):
                            index = bcm.bam_index_path(
                                user, pw, full,
                                gateway_host=gateway_host, target_host=target_host,
                            )
                            if index:
                                bcm.download_to(
                                    user, pw, index,
                                    str(Path(download_dir) / Path(index).name),
                                    gateway_host=gateway_host, target_host=target_host,
                                )
                    register_data_root(download_dir)
                    st.success(
                        f"Downloaded to {target_local}. '{download_dir}' is now "
                        "browsable in Setup."
                    )
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Download failed: {exc}")

        st.caption(
            "Downloading is best for a specific trio or moderate files. Whole-genome "
            "BAMs are large; for those, mount the folder (other tab) and run targeted "
            "analysis in place, or run genome-wide on the server."
        )

        st.divider()
        st.subheader("Fetch a region slice (no server writes, tiny download)")
        st.caption(
            "Reads only the regions you name from a remote BAM (via read-only "
            "`samtools view` on the server) and streams them to a small local BAM "
            "that is indexed and made browsable. Ideal for a targeted gene panel "
            "when you cannot write on the server."
        )
        slice_bam_path = st.text_input(
            "Remote BAM path", st.session_state.get("bcm_slice_bam", ""),
            placeholder="/stornext/.../TrioAnalysis_BH16732/proband.bam",
        )
        gtf_local = st.session_state.get("reference_gtf") or ""
        genes_text = st.text_input(
            "Genes (comma/space separated)", "",
            help="Requires a prepared reference (Setup → download hg38/hg19) so gene "
            "coordinates are known locally. Leave blank to use explicit regions below.",
        )
        regions_text = st.text_input(
            "…or explicit regions", "",
            placeholder="chr3:57192837-57232606 chrX:150000-160000",
        )
        with st.expander("Advanced: samtools location / module"):
            samtools_exe = st.text_input(
                "samtools on server", st.session_state.get("bcm_samtools", "samtools"),
                help="A full path also works, e.g. /opt/conda/envs/bio/bin/samtools.",
            )
            st.session_state.bcm_samtools = samtools_exe
            setup_cmd = st.text_input(
                "Server setup command (optional)", st.session_state.get("bcm_setup", ""),
                placeholder="module load samtools    (or: source ~/.bashrc)",
                help="Runs before samtools. Needed when samtools lives behind a module "
                "or conda env. Everything runs in a login shell so your profile loads.",
            )
            st.session_state.bcm_setup = setup_cmd
            find_cols = st.columns(2)
            if find_cols[0].button("Find samtools on server"):
                try:
                    st.code(bcm.locate_tool(
                        user, pw, samtools_exe.split("/")[-1] or "samtools",
                        gateway_host=gateway_host, target_host=target_host,
                    ))
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Lookup failed: {exc}")
            if find_cols[1].button("Diagnose server (tools + write access)"):
                try:
                    st.code(bcm.diagnose_server(
                        user, pw, gateway_host=gateway_host, target_host=target_host,
                    ))
                except Exception as exc:  # noqa: BLE001
                    st.error(f"Diagnosis failed: {exc}")
        if st.button("Fetch region slice"):
            regions: list[str] = []
            genes = [g.strip() for g in genes_text.replace(",", " ").split() if g.strip()]
            try:
                if genes:
                    if not gtf_local:
                        raise ValueError("Prepare a reference on the Setup page first to resolve gene coordinates.")
                    _, extract, missing = panel_regions(gtf_local, genes, 2000, 5000)
                    regions = extract_to_regions(extract)
                    if missing:
                        st.warning(f"Not found in annotation: {', '.join(missing)}")
                regions += [token for token in regions_text.replace(",", " ").split() if token]
                if not slice_bam_path or not regions:
                    raise ValueError("Provide a remote BAM path and at least one gene or region.")
                out_dir = Path(st.session_state.get("bcm_download_dir", str(Path.home() / "methyl-trio-downloads")))
                local_bam = out_dir / (Path(slice_bam_path).stem + ".slice.bam")
                with st.spinner("Streaming region slice from the server…"):
                    bcm.slice_bam(
                        user, pw, slice_bam_path, regions, str(local_bam),
                        samtools=samtools_exe, setup=st.session_state.get("bcm_setup", ""),
                        gateway_host=gateway_host, target_host=target_host,
                    )
                register_data_root(str(out_dir))
                st.session_state.bcm_slice_bam = slice_bam_path
                st.success(
                    f"Wrote {local_bam} ({len(regions)} regions) and indexed it. "
                    f"'{out_dir}' is now browsable in Setup — select the .slice.bam there."
                )
            except Exception as exc:  # noqa: BLE001
                st.error(f"Slice failed: {exc}")

with tab_mount:
    st.info(
        "This tab never asks for, stores, or transmits your password. You mount the "
        "remote folder yourself in a terminal, where your password and 2FA (Duo) "
        "prompts stay strictly between you and the server. The app only browses the "
        "resulting mount. Each person uses their own account this way."
    )

    with st.expander("One-time setup: install SSHFS"):
        st.markdown(
            "- **macOS:** install [macFUSE](https://osxfuse.github.io/) and "
            "`sshfs` (`brew install --cask macfuse` then `brew install gromgit/fuse/sshfs-mac`), "
            "or use [FUSE-T](https://www.fuse-t.org/) with its sshfs build.\n"
            "- **Linux:** `sudo apt install sshfs` (or your distro's equivalent).\n\n"
            "SSHFS "
            + ("**was detected** on this machine." if remote.sshfs_available()
               else "**was not detected**; install it before mounting.")
        )

    st.subheader("Connection")
    col1, col2 = st.columns(2)
    with col1:
        remote_user = st.text_input("Username", st.session_state.get("remote_user", ""), placeholder="u244415")
        analysis_host = st.text_input(
            "Data host (where the BAMs are)", st.session_state.get("remote_analysis_host", ""),
            placeholder="analysis1.hgsc.bcm.edu",
        )
        remote_path = st.text_input(
            "Remote folder to browse", st.session_state.get("remote_path", ""),
            placeholder="/path/to/your/bams",
            help="The directory on the data host that contains your BAMs.",
        )
    with col2:
        jump_host = st.text_input(
            "Jump / login host (optional)", st.session_state.get("remote_jump_host", ""),
            placeholder="login1.hgsc.bcm.edu",
            help="If you must SSH to a login node first, put it here; the app routes "
            "through it automatically with ProxyJump.",
        )
        jump_user = st.text_input(
            "Jump host username (optional)", st.session_state.get("remote_jump_user", ""),
            placeholder="same as username if blank",
        )
        default_mount = str(Path.home() / "methyl-trio-remote")
        mount_point = st.text_input(
            "Local mount point", st.session_state.get("remote_mount_point", default_mount),
            help="An empty local folder the remote data will appear in.",
        )

    st.session_state.update({
        "remote_user": remote_user, "remote_analysis_host": analysis_host,
        "remote_path": remote_path, "remote_jump_host": jump_host,
        "remote_jump_user": jump_user, "remote_mount_point": mount_point,
    })

    ready = all([remote_user, analysis_host, remote_path, mount_point])
    if ready:
        command = remote.build_sshfs_command(
            remote_user=remote_user, analysis_host=analysis_host, remote_path=remote_path,
            mount_point=mount_point, jump_user=jump_user or None, jump_host=jump_host or None,
            read_only=True,
        )
        st.subheader("Step 1 — mount it in your terminal")
        st.caption(
            "Run this in your own terminal and complete your normal password + Duo "
            "login. It mounts the remote folder read-only."
        )
        st.code(f"mkdir -p {mount_point}\n{remote.command_string(command)}", language="bash")

        st.subheader("Step 2 — use the mount here")
        if st.button("Check mount and use it", type="primary"):
            if remote.is_mounted(mount_point):
                register_data_root(mount_point)
                register_remote_mapping(remote_path, mount_point)
                st.success(
                    f"Mounted. '{mount_point}' is now a browsable location in Setup's "
                    "file pickers. You can paste server paths under "
                    f"'{remote_path}' straight into the picker's 'Go to path' box."
                )
            elif Path(mount_point).is_dir() and any(Path(mount_point).iterdir()):
                register_data_root(mount_point)
                register_remote_mapping(remote_path, mount_point)
                st.warning(
                    "That folder is not detected as a live mount, but it has files, so "
                    "it was added anyway. If browsing fails, re-run the mount command."
                )
            else:
                st.error(
                    "No live mount found at that path yet. Run the command above first, "
                    "finish the password/Duo prompt, then click again."
                )

        with st.expander("Advanced: try to mount from the app (SSH keys/agent only)"):
            st.caption(
                "Only works if you have passwordless SSH keys/agent set up for both "
                "hosts. It will never prompt for a password — if a password is "
                "required it fails quickly and you should use the terminal command."
            )
            if st.button("Attempt key-based mount"):
                if not remote.sshfs_available():
                    st.error("sshfs is not installed on this machine.")
                else:
                    batch_command = remote.build_sshfs_command(
                        remote_user=remote_user, analysis_host=analysis_host,
                        remote_path=remote_path, mount_point=mount_point,
                        jump_user=jump_user or None, jump_host=jump_host or None,
                        read_only=True, batch_mode=True,
                    )
                    ok, output = remote.attempt_mount(batch_command, mount_point)
                    if ok:
                        register_data_root(mount_point)
                        register_remote_mapping(remote_path, mount_point)
                        st.success(f"Mounted and registered '{mount_point}'.")
                    else:
                        st.error("Could not mount without a password. Use the terminal command above.")
                        if output:
                            st.code(output)

        st.subheader("When you're done")
        st.caption("Unmount the remote folder with:")
        st.code(remote.unmount_command(mount_point), language="bash")
    else:
        st.warning("Fill in username, data host, remote folder, and mount point to continue.")

st.divider()
st.caption("Currently browsable locations: " + ", ".join(str(root) for root in data_roots()))
st.info(
    "Analysis reads BAMs directly over the mount, so targeted gene-panel runs are "
    "efficient. Genome-wide runs stream most of each BAM over the network and are "
    "far faster run on the server itself — ask if you'd like a remote-execution "
    "option. Write run outputs to a **local** location on the Run page (the remote "
    "mount is read-only)."
)
