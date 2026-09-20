from __future__ import annotations

import sqlite3
from pathlib import Path

from .config import database_path


SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS imports (
    id TEXT PRIMARY KEY,
    original_filename TEXT NOT NULL,
    archive_sha256 TEXT NOT NULL UNIQUE,
    imported_at TEXT NOT NULL,
    total_members INTEGER NOT NULL,
    imported_records INTEGER NOT NULL,
    skipped_members INTEGER NOT NULL,
    parse_warning_records INTEGER NOT NULL DEFAULT 0,
    parse_warning_count INTEGER NOT NULL DEFAULT 0,
    archive_path TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS sequences (
    id TEXT PRIMARY KEY,
    import_id TEXT NOT NULL REFERENCES imports(id) ON DELETE RESTRICT,
    archive_member_index INTEGER NOT NULL,
    archive_member_name TEXT NOT NULL,
    archive_member_occurrence INTEGER NOT NULL,
    member_sha256 TEXT NOT NULL,
    sequence_sha256 TEXT NOT NULL,
    display_name TEXT NOT NULL,
    length_bp INTEGER NOT NULL,
    topology TEXT NOT NULL CHECK(topology IN ('circular', 'linear', 'unknown')),
    molecule_type TEXT,
    feature_count INTEGER NOT NULL,
    parse_warning_count INTEGER NOT NULL DEFAULT 0,
    parse_warnings_json TEXT NOT NULL DEFAULT '[]',
    raw_genbank TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(import_id, archive_member_index)
);

CREATE INDEX IF NOT EXISTS sequences_display_name_idx ON sequences(display_name);
CREATE INDEX IF NOT EXISTS sequences_sequence_sha_idx ON sequences(sequence_sha256);
"""


def connect(path: Path | None = None) -> sqlite3.Connection:
    db_path = path or database_path()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def initialise(path: Path | None = None) -> None:
    with connect(path) as connection:
        connection.executescript(SCHEMA)
        _add_column_if_missing(
            connection, "imports", "parse_warning_records INTEGER NOT NULL DEFAULT 0"
        )
        _add_column_if_missing(
            connection, "imports", "parse_warning_count INTEGER NOT NULL DEFAULT 0"
        )
        _add_column_if_missing(
            connection, "sequences", "parse_warning_count INTEGER NOT NULL DEFAULT 0"
        )
        _add_column_if_missing(
            connection, "sequences", "parse_warnings_json TEXT NOT NULL DEFAULT '[]'"
        )


def _add_column_if_missing(
    connection: sqlite3.Connection, table: str, definition: str
) -> None:
    """Apply additive migrations for prototype databases created before a new field."""
    column = definition.split()[0]
    columns = {row["name"] for row in connection.execute(f"PRAGMA table_info({table})")}
    if column not in columns:
        connection.execute(f"ALTER TABLE {table} ADD COLUMN {definition}")
