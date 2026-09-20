from __future__ import annotations

import os
from pathlib import Path


def data_dir() -> Path:
    return Path(os.environ.get("MOLBIO_DATA_DIR", "var")).resolve()


def database_path() -> Path:
    return data_dir() / "workbench.sqlite3"


def archive_dir() -> Path:
    return data_dir() / "archives"
