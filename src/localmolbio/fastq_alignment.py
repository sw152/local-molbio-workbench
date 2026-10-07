"""Registered FASTQ -> bounded private snapshots -> local alignment evidence."""
from contextlib import closing
from hashlib import sha256
import gzip
import json
import os
from pathlib import Path
from uuid import UUID, uuid4

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from . import alignment_evidence as core, fastq_validation as registration, job_queue
from .config import data_dir
from .database import connect
from .fastq import FastqError, FastqLimits, inspect_fastq
from .input_validation import InputValidationError, _reference, _hash_managed_file

ADAPTER = 'registered-fastq-minimap2-v1'
SNAPSHOT_BYTES = 32 * 1024**2
LIMITS = FastqLimits(raw_bytes=SNAPSHOT_BYTES, decoded_bytes=SNAPSHOT_BYTES,
                    record_bases=core.MAX_READ, records=core.MAX_READS)


def _tool():
    configured = os.environ.get('MOLBIO_MINIMAP2')
    if not configured:
        raise InputValidationError('minimap2_configuration_missing')
    try:
        return core.tool_identity(configured)
    except (OSError, core.AlignmentError) as exc:
        raise InputValidationError('configured_minimap2_unavailable_or_unsupported') from exc


def _reference_sequence(db, revision):
    row = db.execute('SELECT sequence_text FROM sequence_revisions WHERE id=?', (revision,)).fetchone()
    if row is None:
        raise InputValidationError('reference_missing')
    return row[0]


def _supported(reference, topology, rows):
    try:
        core._validate(reference, ['A'], topology)
    except core.AlignmentError as exc:
        raise InputValidationError(str(exc)) from exc
    if sum(r['size_bytes'] for r in rows) > SNAPSHOT_BYTES:
        raise InputValidationError('alignment_total_original_bytes_limit_exceeded')


def enqueue_alignment(revision_id, input_ids, data_type, idempotency_key, max_attempts=3):
    ids = registration._ids(input_ids)
    if data_type not in core.PRESETS:
        raise InputValidationError('explicit_supported_data_type_required')
    _, tool = _tool()
    with closing(connect()) as db:
        db.execute('BEGIN')
        reference = _reference(db, revision_id)
        rows = registration._rows(db, revision_id, ids)
        identities = [registration._identity(row) for row in rows]
        _supported(_reference_sequence(db, revision_id), reference['topology'], rows)
    return job_queue.enqueue(revision_id, idempotency_key, ADAPTER,
                             {'reference': reference, 'fastq_inputs': identities, 'tool': tool},
                             parameters={'data_type': data_type}, max_attempts=max_attempts)


def execute(claim, pulse):
    try:
        manifest = claim['input_manifest'];expected = manifest['inputs']
        if claim['adapter'] != ADAPTER or manifest['adapter'] != ADAPTER:
            raise InputValidationError('unsupported_adapter_or_parameters')
        if set(claim['parameters']) != {'data_type'} or claim['parameters']['data_type'] not in core.PRESETS:
            raise InputValidationError('unsupported_adapter_or_parameters')
        data_type = claim['parameters']['data_type']
        ids = registration._ids([row['id'] for row in expected['fastq_inputs']])
        identities = {row['id']: row for row in expected['fastq_inputs']}
        if expected['reference'] != manifest['reference']:
            raise InputValidationError('reference_manifest_mismatch')
        revision = manifest['reference']['id']
        job_id = str(UUID(claim['job_id']))
        expected_tool = expected['tool']
    except (KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, InputValidationError):
            raise
        raise InputValidationError('invalid_input_manifest') from exc

    def checked_rows():
        with closing(connect()) as db:
            db.execute('BEGIN')
            if _reference(db, revision) != manifest['reference']:
                raise InputValidationError('reference_changed_since_submission')
            rows = registration._rows(db, revision, ids)
            if any(registration._identity(row) != identities[row['id']] for row in rows):
                raise InputValidationError('fastq_identity_changed_since_submission')
            reference = _reference_sequence(db, revision)
            _supported(reference, manifest['reference']['topology'], rows)
            return reference, rows

    pulse()
    binary, tool = _tool()
    if tool != expected_tool:
        raise InputValidationError('aligner_changed_since_submission')
    reference, rows = checked_rows()
    relative = Path('alignment-attempts') / job_id / str(uuid4())
    attempt = data_dir() / relative
    try:
        (attempt/'inputs').mkdir(parents=True, exist_ok=False)
        (attempt/'attempt.json').write_text(json.dumps({'job_id': job_id, 'attempt_number': claim['attempt'],
            'adapter': ADAPTER, 'reference': manifest['reference'], 'inputs': expected['fastq_inputs'],
            'data_type': data_type, 'tool': expected_tool}, indent=2))
    except OSError as exc:
        raise InputValidationError('alignment_artifact_io_error', retryable=True) from exc
    sequences, sources, snapshots = [], [], []
    decoded_total = bases_total = 0
    try:
        for i, row in enumerate(rows):
            snapshot = attempt/'inputs'/f'{i}.fastq{ ".gz" if row["compression"] == "gzip" else ""}'
            with snapshot.open('xb') as sink:
                _hash_managed_file(row, pulse, data_dir()/'fastq', 'input_outside_managed_fastq',
                                   'input_id', expected_size=row['size_bytes'], sink=sink)
                sink.flush();os.fsync(sink.fileno())
            # Private copy is never overwritten by another attempt; verify raw bytes again.
            if core._file_hash(snapshot) != row['file_sha256']:
                raise InputValidationError('alignment_snapshot_hash_mismatch')
            summary = inspect_fastq(snapshot, row['compression']=='gzip', row['quality_encoding'], LIMITS, pulse=pulse)
            saved = json.loads(row['summary_json'])
            if {k:v for k,v in summary.items() if k!='limits'} != {k:v for k,v in saved.items() if k!='limits'}:
                raise InputValidationError('registered_fastq_summary_does_not_match_bytes')
            decoded_total += summary['decoded_bytes'];bases_total += summary['bases']
            if decoded_total > SNAPSHOT_BYTES or bases_total > core.MAX_BASES or len(sequences)+summary['records'] > core.MAX_READS:
                raise InputValidationError('alignment_aggregate_input_limit_exceeded')
            opener = gzip.open if row['compression']=='gzip' else open
            first_index = len(sequences)
            with opener(snapshot, 'rt', encoding='ascii') as stream:
                for ordinal, (title, sequence, quality) in enumerate(FastqGeneralIterator(stream), 1):
                    if ordinal % 100 == 1:
                        pulse()
                    # Case normalization only. Unsupported IUPAC is rejected by the core.
                    sequence = sequence.upper()
                    sources.append({'query_name': f'q{len(sequences)}', 'input_id': row['id'], 'record_ordinal': ordinal,
                                    'sequence_sha256': sha256(sequence.encode()).hexdigest()})
                    sequences.append(sequence)
            if len(sequences)-first_index != summary['records'] or sum(map(len, sequences[first_index:])) != summary['bases']:
                raise InputValidationError('fastq_parser_count_disagreement')
            if core._file_hash(snapshot) != row['file_sha256']:
                raise InputValidationError('alignment_snapshot_changed_during_parse')
            snapshots.append({'input_id': row['id'], 'relative_path': str(snapshot.relative_to(attempt)),
                              'sha256': row['file_sha256']})
        core._validate(reference, sequences, manifest['reference']['topology'])
        source_text = json.dumps(sources, indent=2)
        source_hash = sha256(source_text.encode()).hexdigest()
        (attempt/'sources.json').write_text(source_text)
        result = core.align(reference, sequences, topology=manifest['reference']['topology'], data_type=data_type,
                            binary=binary, output_dir=attempt/'alignment', pulse=pulse)
        if result['tool'] != expected_tool:
            raise InputValidationError('aligner_changed_since_submission')
        # Conservative publication policy: original registered bytes must still match too.
        for row, snapshot in zip(rows, snapshots):
            pulse()
            _hash_managed_file(row, pulse, data_dir()/'fastq', 'input_outside_managed_fastq',
                               'input_id', expected_size=row['size_bytes'])
            if core._file_hash(attempt/snapshot['relative_path']) != snapshot['sha256']:
                raise InputValidationError('alignment_snapshot_changed_during_run')
        if core._file_hash(attempt/'sources.json') != source_hash:
            raise InputValidationError('alignment_source_mapping_changed')
        if checked_rows() != (reference, rows):
            raise InputValidationError('registered_inputs_changed_during_alignment')
        result = {**result, 'schema': 'localmolbio.registered-fastq-alignment',
                  'scope': 'registered_fastq_local_alignments_only', 'analysis_performed': True,
                  'job_id': job_id, 'attempt_number': claim['attempt'],
                  'reference': manifest['reference'], 'input_identities': expected['fastq_inputs'],
                  'artifact_directory': str(relative), 'snapshots': snapshots,
                  'sources_sha256': source_hash, 'sources': sources}
        pulse()
        # Diagnostic attempt result only. The queue's finish token gates publication.
        (attempt/'attempt-result.json').write_text(json.dumps(result, indent=2))
        return result
    except (FastqError, core.AlignmentError) as exc:
        raise InputValidationError(str(exc)) from exc
    except OSError as exc:
        raise InputValidationError('alignment_artifact_io_error', retryable=True) from exc
