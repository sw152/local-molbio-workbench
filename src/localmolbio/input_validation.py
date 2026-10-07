"""Point-in-time identity checks for already registered AB1 inputs; no alignment."""
from contextlib import closing
from hashlib import sha256
import os
from pathlib import Path
import stat

from .config import read_dir
from .database import connect
from . import job_queue

ADAPTER = 'registered-ab1-input-check-v1'


class InputValidationError(ValueError):
    def __init__(self, code, retryable=False):
        super().__init__(code)
        self.retryable = retryable


def _reference(db, revision):
    row = db.execute('SELECT id,sequence_text,sequence_sha256,topology FROM sequence_revisions WHERE id=?', (revision,)).fetchone()
    if row is None:
        raise InputValidationError('reference_missing')
    if sha256(row['sequence_text'].encode()).hexdigest() != row['sequence_sha256']:
        raise InputValidationError('reference_hash_mismatch')
    return {key: row[key] for key in ('id', 'sequence_sha256', 'topology')}


def _reads(db, revision, ids):
    rows = []
    for read_id in ids:
        row = db.execute('SELECT id,sequence_revision_id,file_sha256,storage_path,file_format FROM sequencing_reads WHERE id=?', (read_id,)).fetchone()
        if row is None or row['sequence_revision_id'] != revision:
            raise InputValidationError('read_missing_or_wrong_revision')
        if row['file_format'] != 'ab1':
            raise InputValidationError('unsupported_input_format')
        rows.append(dict(row))
    return rows


def _ids(ids):
    if not isinstance(ids, list) or not 1 <= len(ids) <= 100 or any(not isinstance(i, str) or not i for i in ids) or len(set(ids)) != len(ids):
        raise InputValidationError('expected_1_to_100_distinct_read_ids')
    return sorted(ids)


def enqueue_check(revision_id, read_ids, idempotency_key, max_attempts=3):
    ids = _ids(read_ids)
    with closing(connect()) as db:
        db.execute('BEGIN')
        reference = _reference(db, revision_id)
        rows = _reads(db, revision_id, ids)
    # Include the checked reference in the immutable inputs too. If it changes between
    # this snapshot and enqueue, the adapter will reject the inconsistent request.
    inputs = {'reference': reference, 'reads': [{'id': r['id'], 'sha256': r['file_sha256']} for r in rows]}
    return job_queue.enqueue(revision_id, idempotency_key, ADAPTER, inputs, max_attempts=max_attempts)


def _signature(info):
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def _hash_managed_file(row, pulse, directory, outside_code, identity_key, expected_size=None, sink=None):
    path = Path(row['storage_path'])
    try:
        if path.is_symlink():
            raise InputValidationError('symlink_input_rejected')
        resolved = path.resolve(strict=True)
        if not resolved.is_relative_to(directory.resolve()):
            raise InputValidationError(outside_code)
        # NONBLOCK prevents an unexpectedly replaced FIFO from blocking before fstat.
        fd = os.open(resolved, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        try:
            before = os.fstat(fd)
            if not stat.S_ISREG(before.st_mode):
                raise InputValidationError('input_not_regular_file')
            if expected_size is not None and before.st_size != expected_size:
                raise InputValidationError('input_size_mismatch')
            digest = sha256()
            copied_bytes = 0
            with os.fdopen(fd, 'rb', closefd=False) as stream:
                while True:
                    pulse()
                    block = stream.read(1024 * 1024)
                    if not block:
                        break
                    copied_bytes += len(block)
                    if expected_size is not None and copied_bytes > expected_size:
                        raise InputValidationError('input_changed_during_check')
                    digest.update(block)
                    if sink is not None:
                        sink.write(block)
            after = os.fstat(fd)
            if path.is_symlink() or _signature(before) != _signature(after) or _signature(after) != _signature(path.stat()):
                raise InputValidationError('input_changed_during_check')
        finally:
            os.close(fd)
        if digest.hexdigest() != row['file_sha256']:
            raise InputValidationError('input_hash_mismatch')
        return {identity_key: row['id'], 'sha256': digest.hexdigest(), 'size_bytes': after.st_size}
    except FileNotFoundError as exc:
        raise InputValidationError('input_file_missing') from exc
    except PermissionError as exc:
        raise InputValidationError('input_not_readable') from exc
    except OSError as exc:
        raise InputValidationError('input_io_error', retryable=True) from exc


def _hash_registered_file(row, pulse):
    return _hash_managed_file(row, pulse, read_dir(), 'input_outside_managed_reads', 'read_id')


def validate_inputs(claim, pulse):
    manifest = claim['input_manifest']
    try:
        if claim['adapter'] != ADAPTER or manifest['adapter'] != ADAPTER or claim['parameters'] != {}:
            raise InputValidationError('unsupported_adapter_or_parameters')
        expected = manifest['inputs']
        entries = expected['reads']
        ids = _ids([r['id'] for r in entries])
        hashes = {r['id']: r['sha256'] for r in entries}
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
            rows = _reads(db, revision, ids)
            if any(r['file_sha256'] != hashes[r['id']] for r in rows):
                raise InputValidationError('read_identity_changed_since_submission')
            return rows

    pulse()
    rows = checked_rows()
    files = [_hash_registered_file(row, pulse) for row in rows]
    pulse()
    if checked_rows() != rows:
        raise InputValidationError('registered_inputs_changed_during_check')
    return {'schema': 'localmolbio.input-check', 'schema_version': 1, 'scope': 'registered_input_integrity_only',
            'adapter': ADAPTER, 'reference': manifest['reference'], 'files': files,
            'analysis_performed': False, 'whole_reference_verified': False}
