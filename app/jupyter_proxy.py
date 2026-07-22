"""Expose MANGO through jupyter-server-proxy inside a Terra Cloud Environment.

AnVIL / Terra interactive environments are launched by Leonardo as a Jupyter
server. ``jupyter-server-proxy`` lets a Jupyter server front arbitrary web apps;
registering this ``launch`` callable under the ``jupyter_serverproxy_servers``
entry point (see ``pyproject.toml``) makes MANGO appear in the Jupyter Launcher
and be served at ``.../proxy/mango/`` — so the same Streamlit UI runs entirely
inside the AnVIL security and billing boundary.
"""
from __future__ import annotations

from pathlib import Path


def _app_path() -> str:
    return str(Path(__file__).with_name("streamlit_app.py"))


def launch() -> dict:
    """Return the jupyter-server-proxy configuration that starts Streamlit.

    ``{port}`` and ``{base_url}`` are substituted by jupyter-server-proxy at
    launch. Streamlit is told its base path so assets and the websocket resolve
    correctly behind the Jupyter proxy.
    """
    return {
        "command": [
            "streamlit",
            "run",
            _app_path(),
            "--server.port={port}",
            "--server.address=127.0.0.1",
            "--server.headless=true",
            "--server.enableCORS=false",
            "--server.enableXsrfProtection=false",
            "--server.baseUrlPath={base_url}proxy/mango",
        ],
        "timeout": 120,
        "absolute_url": False,
        "launcher_entry": {
            "title": "MANGO",
            "enabled": True,
        },
    }
