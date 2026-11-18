"""Launch the ASE database web viewer with a local workaround for metadata.

Usage:
    source /home/caterine/psi4conda/bin/activate dft_env
    cd /home/caterine/Documentos/tesis/DFT
    python abrir_db_firefox.py
"""

from __future__ import annotations

import argparse
import webbrowser

from ase.db import connect
import ase.db.app as app


DEFAULT_DB = "bio_molecules_smiles_estereoquimica.db"


def launch_database(db_path: str) -> None:
    db = connect(db_path)
    if getattr(db, "_metadata", None) is None:
        db._metadata = {}

    db_app = app.DBApp()
    db_app.add_project("default", db)

    try:
        webbrowser.get("firefox").open_new_tab("http://127.0.0.1:5000")
    except Exception:
        webbrowser.open_new_tab("http://127.0.0.1:5000")

    db_app.flask.run(host="0.0.0.0", debug=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Open an ASE database in the browser."
    )
    parser.add_argument(
        "db_path",
        nargs="?",
        default=DEFAULT_DB,
        help=f"Database file to open (default: {DEFAULT_DB}).",
    )
    args = parser.parse_args()
    launch_database(args.db_path)


if __name__ == "__main__":
    main()