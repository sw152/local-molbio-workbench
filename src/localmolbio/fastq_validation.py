"""Identity-only queue adapter for registered FASTQ bytes and declared metadata."""
from contextlib import closing
from hashlib import sha256
import json
import re

from .config import data_dir
from .database import connect
from .fastq import FastqLimits
from .input_validation import InputValidationError, _reference, _hash_managed_file
from . import job_queue

ADAPTER = 'registered-fastq-input-check-v1'


def _ids(ids):
    if not isinstance(ids, list) or not 1 <= len(ids) <= 100 or any(not isinstance(i, str) or not i for i in ids) or len(set(ids)) != len(ids):
        raise InputValidationError('expected_1_to_100_distinct_fastq_ids')
    return sorted(ids)


def _identity(row):
    # Bind the saved summary too, but do not claim that it was recomputed by this check.
    try:
        summary = json.loads(row['summary_json'])
        if (not isinstance(summary, dict) or summary.get('schema_version') != 1
                or summary.get('parser') != 'bounded-fastq-v1'
                or summary.get('scope') != 'fastq_format_and_quality_summary'
                or summary.get('compression') != row['compression']
                or summary.get('quality_encoding') != row['quality_encoding']
                or summary.get('quality_encoding_source') != 'explicit_user_declaration'
                or summary.get('read_layout') != 'single_file_unpaired'
                or summary.get('analysis_performed') is not False
                or summary.get('whole_reference_verified') is not False):
            raise InputValidationError('unsupported_fastq_registration_metadata')
        canonical = json.dumps(summary, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)
    except (ValueError, TypeError) as exc:
        if isinstance(exc, InputValidationError):
            raise
        raise InputValidationError('invalid_fastq_registration_summary') from exc
    if (row['compression'] not in ('none', 'gzip') or row['quality_encoding'] != 'phred33'
            or not isinstance(row['file_sha256'], str) or not re.fullmatch('[0-9a-f]{64}', row['file_sha256'])
            or type(row['size_bytes']) is not int or not 0 < row['size_bytes'] <= FastqLimits().raw_bytes):
        raise InputValidationError('unsupported_fastq_registration_metadata')
    return {'id': row['id'], 'sha256': row['file_sha256'], 'size_bytes': row['size_bytes'],
            'compression': row['compression'], 'quality_encoding': row['quality_encoding'],
            'summary_sha256': sha256(canonical.encode()).hexdigest()}


def _rows(db, revision, ids):
    rows = []
    for identity in ids:
        row = db.execute('SELECT * FROM fastq_inputs WHERE id=?', (identity,)).fetchone()
        if row is None or row['sequence_revision_id'] != revision:
            raise InputValidationError('fastq_missing_or_wrong_revision')
        rows.append(dict(row))
    return rows


def enqueue_check(revision_id, input_ids, idempotency_key, max_attempts=3):
    ids = _ids(input_ids)
    with closing(connect()) as db:
        db.execute('BEGIN')
        reference = _reference(db, revision_id)
        inputs = [_identity(row) for row in _rows(db, revision_id, ids)]
    return job_queue.enqueue(revision_id, idempotency_key, ADAPTER,
                             {'reference': reference, 'fastq_inputs': inputs}, max_attempts=max_attempts)


def validate_inputs(claim, pulse):
    try:
        manifest = claim['input_manifest']
        if claim['adapter'] != ADAPTER or manifest['adapter'] != ADAPTER or claim['parameters'] != {}:
            raise InputValidationError('unsupported_adapter_or_parameters')
        expected = manifest['inputs']
        entries = expected['fastq_inputs']
        ids = _ids([r['id'] for r in entries])
        identities = {r['id']: r for r in entries}
        if expected['reference'] != manifest['reference']:
            raise InputValidationError('reference_manifest_mismatch')
        revision = manifest['reference']['id']
    except (KeyError, TypeError) as exc:
        raise InputValidationError('invalid_input_manifest') from exc

    def checked_rows():
        with closing(connect()) as db:
            db.execute('BEGIN')
            if _reference(db, revision) != manifest['reference']:
                raise InputValidationError('reference_changed_since_submission')
            rows = _rows(db, revision, ids)
            if any(_identity(row) != identities[row['id']] for row in rows):
                raise InputValidationError('fastq_identity_changed_since_submission')
            return rows

    pulse()
    rows = checked_rows()
    files = []
    for row in rows:
        checked = _hash_managed_file(row, pulse, data_dir() / 'fastq', 'input_outside_managed_fastq',
                                     'input_id', expected_size=row['size_bytes'])
        files.append({**checked, **{key: identities[row['id']][key]
                     for key in ('compression', 'quality_encoding', 'summary_sha256')}})
    pulse()
    if checked_rows() != rows:
        raise InputValidationError('registered_inputs_changed_during_check')
    return {'schema': 'localmolbio.fastq-input-check', 'schema_version': 1,
            'scope': 'registered_fastq_integrity_only', 'adapter': ADAPTER,
            'reference': manifest['reference'], 'files': files, 'summary_recomputed': False,
            'analysis_performed': False, 'whole_reference_verified': False}
