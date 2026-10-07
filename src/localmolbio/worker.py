"""Explicit one-shot local worker: python -m localmolbio.worker --help."""
import argparse
import json
import socket
from uuid import uuid4

from .database import initialise
from . import job_queue, fastq_validation, fastq_alignment
from .input_validation import ADAPTER, InputValidationError, enqueue_check, validate_inputs


def run_once(worker_id, lease_seconds=60, adapter=ADAPTER):
    if adapter not in (ADAPTER, fastq_validation.ADAPTER, fastq_alignment.ADAPTER):
        raise ValueError('Unsupported worker adapter')
    validator = {ADAPTER: validate_inputs, fastq_validation.ADAPTER: fastq_validation.validate_inputs,
                 fastq_alignment.ADAPTER: fastq_alignment.execute}[adapter]
    claim = job_queue.claim(worker_id, adapter, lease_seconds)
    if claim is None:
        return {'status': 'idle'}
    job, token = claim['job_id'], claim['lease_token']
    response = {'job_id': job, 'attempt': claim['attempt']}
    try:
        result = validator(claim, lambda: job_queue.renew(job, token, lease_seconds))
        job_queue.finish(job, token, result)
        return {**response, 'status': 'succeeded', 'scope': result['scope'], 'analysis_performed': result.get('analysis_performed', False)}
    except job_queue.QueueConflict:
        return {**response, 'status': 'lease_lost'}
    except Exception as exc:
        # Never persist arbitrary exception text: filesystem paths and input contents
        # can be sensitive. Unknown failures are terminal and preserve attempt identity.
        code = str(exc) if isinstance(exc, InputValidationError) else ('unexpected_alignment_error' if adapter == fastq_alignment.ADAPTER else 'unexpected_input_check_error')
        retryable = isinstance(exc, InputValidationError) and exc.retryable
        try:
            job_queue.fail(job, token, code, retryable=retryable)
        except job_queue.QueueConflict:
            return {**response, 'status': 'lease_lost'}
        return {**response, 'status': 'attempt_failed', 'error': code}


def main():
    parser = argparse.ArgumentParser(description='Registered input checks and explicitly selected bounded FASTQ local alignments; no consensus or plasmid pass verdict.')
    sub = parser.add_subparsers(dest='command', required=True)
    enqueue = sub.add_parser('enqueue-check', help='Queue a check of previously uploaded AB1 read IDs')
    enqueue.add_argument('--revision', required=True)
    enqueue.add_argument('--read-id', action='append', required=True)
    enqueue.add_argument('--key', required=True, help='Stable unique submission key')
    enqueue.add_argument('--max-attempts', type=int, choices=range(1, 11), default=3)
    fastq = sub.add_parser('enqueue-fastq-check', help='Queue an identity check of registered FASTQ input IDs')
    fastq.add_argument('--revision', required=True)
    fastq.add_argument('--input-id', action='append', required=True)
    fastq.add_argument('--key', required=True)
    fastq.add_argument('--max-attempts', type=int, choices=range(1, 11), default=3)
    alignment = sub.add_parser('enqueue-alignment', help='Queue bounded local alignment of registered FASTQ inputs')
    alignment.add_argument('--revision', required=True)
    alignment.add_argument('--input-id', action='append', required=True)
    alignment.add_argument('--data-type', choices=tuple(fastq_alignment.core.PRESETS), required=True)
    alignment.add_argument('--key', required=True)
    alignment.add_argument('--max-attempts', type=int, choices=range(1, 11), default=3)
    once = sub.add_parser('once', help='Claim at most one supported task, then exit')
    once.add_argument('--worker-id', default=f'{socket.gethostname()}-{uuid4()}')
    once.add_argument('--adapter', choices=('ab1', 'fastq', 'alignment'), default='ab1', help='Only claim this adapter (default: ab1)')
    once.add_argument('--lease-seconds', type=int, default=60)
    args = parser.parse_args()
    initialise()
    try:
        if args.command == 'enqueue-check':
            result = {'job_id': enqueue_check(args.revision, args.read_id, args.key, args.max_attempts), 'status': 'submitted'}
        elif args.command == 'enqueue-fastq-check':
            result = {'job_id': fastq_validation.enqueue_check(args.revision, args.input_id, args.key, args.max_attempts), 'status': 'submitted'}
        elif args.command == 'enqueue-alignment':
            result = {'job_id': fastq_alignment.enqueue_alignment(args.revision, args.input_id, args.data_type, args.key, args.max_attempts), 'status': 'submitted'}
        else:
            adapter = {'ab1': ADAPTER, 'fastq': fastq_validation.ADAPTER, 'alignment': fastq_alignment.ADAPTER}[args.adapter]
            result = run_once(args.worker_id, args.lease_seconds, adapter)
    except (ValueError, job_queue.QueueConflict) as exc:
        print(json.dumps({'status': 'rejected', 'error': str(exc)}))
        return 2
    print(json.dumps(result))
    return 0 if result['status'] in {'submitted', 'idle', 'succeeded'} else 1


if __name__ == '__main__':
    raise SystemExit(main())
