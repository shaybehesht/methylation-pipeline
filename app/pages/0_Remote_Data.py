from pathlib import Path

import streamlit as st

from app import remote
from app.file_picker import data_roots, register_data_root
from app.state import initialize

initialize()
st.title("0. Remote data (SSH / SSHFS)")
st.caption(
    "Browse BAMs that live on a remote server by mounting the remote folder over "
    "SSH. Files are read on demand — nothing is permanently downloaded."
)

st.info(
    "This app never asks for, stores, or transmits your password. You mount the "
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
            st.success(
                f"Mounted. '{mount_point}' is now a browsable location in Setup's "
                "file pickers."
            )
        elif Path(mount_point).is_dir() and any(Path(mount_point).iterdir()):
            register_data_root(mount_point)
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
