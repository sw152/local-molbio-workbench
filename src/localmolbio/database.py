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
    features_json TEXT NOT NULL DEFAULT '[]',
    parse_warning_count INTEGER NOT NULL DEFAULT 0,
    parse_warnings_json TEXT NOT NULL DEFAULT '[]',
    raw_genbank TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(import_id, archive_member_index)
);

CREATE INDEX IF NOT EXISTS sequences_display_name_idx ON sequences(display_name);
CREATE INDEX IF NOT EXISTS sequences_sequence_sha_idx ON sequences(sequence_sha256);

CREATE TABLE IF NOT EXISTS sequence_revisions (
    id TEXT PRIMARY KEY,
    sequence_id TEXT NOT NULL REFERENCES sequences(id) ON DELETE RESTRICT,
    parent_revision_id TEXT REFERENCES sequence_revisions(id) ON DELETE RESTRICT,
    revision_number INTEGER NOT NULL,
    label TEXT NOT NULL,
    sequence_text TEXT NOT NULL,
    sequence_sha256 TEXT NOT NULL,
    topology TEXT NOT NULL CHECK(topology IN ('circular', 'linear', 'unknown')),
    features_json TEXT NOT NULL DEFAULT '[]',
    source_kind TEXT NOT NULL CHECK(source_kind IN ('import', 'manual-edit')),
    created_at TEXT NOT NULL,
    UNIQUE(sequence_id, revision_number)
);

CREATE INDEX IF NOT EXISTS sequence_revisions_sequence_idx ON sequence_revisions(sequence_id, revision_number DESC);

CREATE TABLE IF NOT EXISTS primers (
    id TEXT PRIMARY KEY,
    sequence_revision_id TEXT NOT NULL REFERENCES sequence_revisions(id) ON DELETE RESTRICT,
    name TEXT,
    sequence_text TEXT NOT NULL,
    direction TEXT NOT NULL CHECK(direction IN ('forward', 'reverse')),
    purpose TEXT NOT NULL,
    binding_start INTEGER,
    binding_end INTEGER,
    metrics_json TEXT NOT NULL DEFAULT '{}',
    design_parameters_json TEXT NOT NULL DEFAULT '{}',
    selection_state TEXT NOT NULL CHECK(selection_state IN ('candidate', 'selected', 'archived')),
    created_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS primers_revision_idx ON primers(sequence_revision_id, selection_state);

CREATE TABLE IF NOT EXISTS analysis_jobs (
    id TEXT PRIMARY KEY,
    sequence_revision_id TEXT REFERENCES sequence_revisions(id) ON DELETE RESTRICT,
    job_kind TEXT NOT NULL CHECK(job_kind IN ('primer-design', 'sanger-verification', 'plasmid-verification')),
    status TEXT NOT NULL CHECK(status IN ('queued', 'running', 'succeeded', 'failed', 'cancelled')),
    parameters_json TEXT NOT NULL DEFAULT '{}',
    input_manifest_json TEXT NOT NULL DEFAULT '{}',
    result_summary_json TEXT NOT NULL DEFAULT '{}',
    error_detail TEXT,
    created_at TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT
);

CREATE INDEX IF NOT EXISTS analysis_jobs_status_idx ON analysis_jobs(status, created_at DESC);

CREATE TABLE IF NOT EXISTS sequencing_reads (
    id TEXT PRIMARY KEY,
    sequence_revision_id TEXT NOT NULL REFERENCES sequence_revisions(id) ON DELETE RESTRICT,
    original_filename TEXT NOT NULL,
    file_sha256 TEXT NOT NULL,
    storage_path TEXT NOT NULL,
    file_format TEXT NOT NULL CHECK(file_format IN ('ab1', 'fastq')),
    direction TEXT NOT NULL CHECK(direction IN ('forward', 'reverse', 'unknown')),
    base_sequence TEXT NOT NULL,
    length_bp INTEGER NOT NULL,
    quality_summary_json TEXT NOT NULL DEFAULT '{}',
    parser_metadata_json TEXT NOT NULL DEFAULT '{}',
    created_at TEXT NOT NULL,
    UNIQUE(sequence_revision_id, file_sha256)
);

CREATE INDEX IF NOT EXISTS sequencing_reads_revision_idx ON sequencing_reads(sequence_revision_id, created_at DESC);

CREATE TABLE IF NOT EXISTS audit_events (
    id TEXT PRIMARY KEY,
    object_type TEXT NOT NULL,
    object_id TEXT NOT NULL,
    action TEXT NOT NULL,
    payload_json TEXT NOT NULL DEFAULT '{}',
    created_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS audit_events_object_idx ON audit_events(object_type, object_id, created_at DESC);
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
            connection, "sequences", "features_json TEXT NOT NULL DEFAULT '[]'"
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
