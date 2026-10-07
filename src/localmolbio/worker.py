"""Explicit one-shot local worker: python -m localmolbio.worker --help."""
import argparse
import json
import socket
from uuid import uuid4

from .database import initialise
from . import job_queue
from .input_validation import ADAPTER, InputValidationError, enqueue_check, validate_inputs


def run_once(worker_id, lease_seconds=60):
    claim = job_queue.claim(worker_id, ADAPTER, lease_seconds)
    if claim is None:
        return {'status': 'idle'}
    job, token = claim['job_id'], claim['lease_token']
    response = {'job_id': job, 'attempt': claim['attempt']}
    try:
        result = validate_inputs(claim, lambda: job_queue.renew(job, token, lease_seconds))
        job_queue.finish(job, token, result)
        return {**response, 'status': 'succeeded', 'scope': result['scope'], 'analysis_performed': False}
    except job_queue.QueueConflict:
        return {**response, 'status': 'lease_lost'}
    except Exception as exc:
        # Never persist arbitrary exception text: filesystem paths and input contents
        # can be sensitive. Unknown failures are terminal and preserve attempt identity.
        code = str(exc) if isinstance(exc, InputValidationError) else 'unexpected_input_check_error'
        retryable = isinstance(exc, InputValidationError) and exc.retryable
        try:
            job_queue.fail(job, token, code, retryable=retryable)
        except job_queue.QueueConflict:
            return {**response, 'status': 'lease_lost'}
        return {**response, 'status': 'attempt_failed', 'error': code}


def main():
    parser = argparse.ArgumentParser(description='Registered AB1 identity checks only; does not perform sequencing analysis.')
    sub = parser.add_subparsers(dest='command', required=True)
    enqueue = sub.add_parser('enqueue-check', help='Queue a check of previously uploaded AB1 read IDs')
    enqueue.add_argument('--revision', required=True)
    enqueue.add_argument('--read-id', action='append', required=True)
    enqueue.add_argument('--key', required=True, help='Stable unique submission key')
    enqueue.add_argument('--max-attempts', type=int, choices=range(1, 11), default=3)
    once = sub.add_parser('once', help='Claim at most one supported task, then exit')
    once.add_argument('--worker-id', default=f'{socket.gethostname()}-{uuid4()}')
    once.add_argument('--lease-seconds', type=int, default=60)
    args = parser.parse_args()
    initialise()
    try:
        if args.command == 'enqueue-check':
            result = {'job_id': enqueue_check(args.revision, args.read_id, args.key, args.max_attempts), 'status': 'submitted'}
        else:
            result = run_once(args.worker_id, args.lease_seconds)
    except (ValueError, job_queue.QueueConflict) as exc:
        print(json.dumps({'status': 'rejected', 'error': str(exc)}))
        return 2
    print(json.dumps(result))
    return 0 if result['status'] in {'submitted', 'idle', 'succeeded'} else 1


if __name__ == '__main__':
    raise SystemExit(main())
