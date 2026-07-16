"""Console entry point so the app starts with a single command.

After ``pip install -e .`` the ``methyl-trio`` command launches Streamlit with
the correct app path, regardless of the current working directory.
"""
from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    from streamlit.web import cli as stcli

    app_path = Path(__file__).with_name("streamlit_app.py")
    sys.argv = ["streamlit", "run", str(app_path), *sys.argv[1:]]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
