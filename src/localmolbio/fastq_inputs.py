"""Atomic registration of original FASTQ bytes, separate from Sanger read records."""
from contextlib import closing
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4

from .config import data_dir
from .database import connect
from .fastq import FastqError, FastqLimits, inspect_fastq

SCHEMA = """
CREATE TABLE IF NOT EXISTS fastq_inputs (
    id TEXT PRIMARY KEY,
    sequence_revision_id TEXT NOT NULL REFERENCES sequence_revisions(id) ON DELETE RESTRICT,
    original_filename TEXT NOT NULL,
    file_sha256 TEXT NOT NULL,
    storage_path TEXT NOT NULL,
    size_bytes INTEGER NOT NULL,
    compression TEXT NOT NULL CHECK(compression IN ('none','gzip')),
    quality_encoding TEXT NOT NULL CHECK(quality_encoding='phred33'),
    summary_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(sequence_revision_id,file_sha256)
);
CREATE INDEX IF NOT EXISTS fastq_revision_idx ON fastq_inputs(sequence_revision_id,created_at DESC);
"""


class FastqConflict(FastqError):
    pass


def _hash(path):
    digest = sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(1024**2):
            digest.update(chunk)
    return digest.hexdigest()


def register_fastq(revision, filename, source, quality_encoding, limits=FastqLimits()):
    if quality_encoding != 'phred33':
        raise FastqError('quality_encoding_must_be_explicit_phred33')
    if not isinstance(filename, str) or not filename or len(filename) > 255:
        raise FastqError('invalid_filename')
    name = filename.lower()
    if not name.endswith(('.fastq','.fq','.fastq.gz','.fq.gz')):
        raise FastqError('expected_fastq_or_fastq_gz')
    compressed = name.endswith('.gz')
    with closing(connect()) as db:
        if db.execute('SELECT 1 FROM sequence_revisions WHERE id=?',(revision,)).fetchone() is None:
            raise FastqError('reference_missing')
    directory = data_dir() / 'fastq'
    staging = directory / '.incoming'
    staging.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(dir=staging, prefix='upload-', delete=False) as output:
            temporary = Path(output.name)
            digest = sha256(); size = 0
            while chunk := source.read(1024**2):
                size += len(chunk)
                if size > limits.raw_bytes:
                    raise FastqError('raw_size_limit_exceeded')
                digest.update(chunk); output.write(chunk)
            output.flush(); os.fsync(output.fileno())
        summary = inspect_fastq(temporary, compressed, quality_encoding, limits)
        file_hash = digest.hexdigest()
        destination = directory / (file_hash + ('.fastq.gz' if compressed else '.fastq'))
        now = datetime.now(timezone.utc).isoformat(); identity = str(uuid4())
        with closing(connect()) as db:
            try:
                db.execute('BEGIN IMMEDIATE')
                if db.execute('SELECT 1 FROM sequence_revisions WHERE id=?',(revision,)).fetchone() is None:
                    raise FastqError('reference_missing')
                if db.execute('SELECT 1 FROM fastq_inputs WHERE sequence_revision_id=? AND file_sha256=?',(revision,file_hash)).fetchone():
                    raise FastqConflict('fastq_already_registered_for_revision')
                try:
                    os.link(temporary, destination)
                except FileExistsError:
                    if destination.is_symlink() or not destination.is_file() or _hash(destination) != file_hash:
                        raise FastqConflict('existing_fastq_storage_integrity_mismatch')
                db.execute('INSERT INTO fastq_inputs VALUES (?,?,?,?,?,?,?,?,?,?)',
                           (identity,revision,filename,file_hash,str(destination),size,summary['compression'],quality_encoding,json.dumps(summary),now))
                db.execute('INSERT INTO audit_events VALUES (?,?,?,?,?,?)',
                           (str(uuid4()),'fastq-input',identity,'registered',json.dumps({'sequence_revision_id':revision,'file_sha256':file_hash,'size_bytes':size}),now))
                db.commit()
            except BaseException:
                db.rollback()
                # Keep any content-addressed file: another revision may reference it.
                raise
        return {'id':identity,'sequence_revision_id':revision,'original_filename':filename,
                'file_sha256':file_hash,'size_bytes':size,'compression':summary['compression'],
                'quality_encoding':quality_encoding,'summary':summary,'created_at':now}
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
