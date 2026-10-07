"""Revision-scoped alignment jobs; publish only queue-committed evidence."""
from contextlib import closing
import json
from typing import Literal, Optional

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, ConfigDict, Field

from .database import connect
from .fastq_alignment import ADAPTER, enqueue_alignment
from .input_validation import InputValidationError
from .job_queue import QueueConflict, cancel

router = APIRouter(prefix='/api/revisions/{revision_id}/alignments')
SCOPE = 'registered_fastq_local_alignments_only'
# Explicit public projections: paths, raw sequences and lease credentials stay internal.
RESULT_FIELDS = ('schema', 'schema_version', 'scope', 'analysis_performed', 'attempt_number',
                 'tool', 'data_type', 'topology', 'coordinate_convention', 'reference_length_bp',
                 'reference_sha256', 'limits', 'secondary_limit', 'candidate_search_exhaustive',
                 'base_quality_used', 'pairing_used', 'circular_mapq_recalibrated',
                 'consensus_performed', 'whole_reference_verified', 'sources_sha256')
READ_FIELDS = ('read_index', 'length_bp', 'alignments', 'withheld', 'status',
               'uniqueness_assessed', 'requires_review')


class Submission(BaseModel):
    model_config = ConfigDict(extra='forbid')
    input_ids: list[str] = Field(min_length=1, max_length=100)
    data_type: Literal['ont-noisy', 'ont-high-accuracy', 'pacbio-hifi', 'short-single']
    idempotency_key: str = Field(min_length=1, max_length=200)
    max_attempts: int = Field(default=3, ge=1, le=10, strict=True)


def _revision(db, revision):
    if db.execute('SELECT 1 FROM sequence_revisions WHERE id=?', (revision,)).fetchone() is None:
        raise HTTPException(404, 'Sequence revision not found')


def _job(db, revision, job):
    row = db.execute('''SELECT j.id,j.status,j.created_at,j.started_at,j.completed_at,j.error_detail,
                      j.input_manifest_json,j.parameters_json,j.result_summary_json,
                      q.attempt_count,q.max_attempts
                      FROM analysis_jobs j JOIN queued_jobs q ON q.job_id=j.id
                      WHERE j.id=? AND j.sequence_revision_id=? AND q.adapter=?''',
                     (job, revision, ADAPTER)).fetchone()
    if row is None:
        raise HTTPException(404, 'Alignment not found for this revision')
    return row


def _summary(row):
    result = {k:row[k] for k in ('id','status','created_at','started_at','completed_at',
                               'error_detail','attempt_count','max_attempts')}
    manifest = json.loads(row['input_manifest_json'])
    result.update(input_count=len(manifest['inputs']['fastq_inputs']),
                  data_type=json.loads(row['parameters_json'])['data_type'], scope=SCOPE)
    return result


def _page(items, total, limit, offset):
    return dict(items=items, total=total, limit=limit, offset=offset, has_more=offset+len(items)<total)


def _published(row):
    if row['status'] != 'succeeded':
        raise HTTPException(409, 'alignment_evidence_not_published')
    result = json.loads(row['result_summary_json'])
    if result.get('schema') != 'localmolbio.registered-fastq-alignment' or result.get('schema_version') != 1:
        raise HTTPException(409, 'unsupported_alignment_evidence_schema')
    return result


def _validate_sources(evidence):
    # Fail closed on broken provenance, rather than displaying another record's evidence.
    try:
        all_reads, sources = evidence['reads'], evidence['sources']
        identities = {i['id'] for i in evidence['input_identities']}
        ordinals = {i:0 for i in identities}
        valid = len(all_reads) == len(sources) == len(evidence['query_sha256'])
        # sources_sha256 hashes the original file bytes, not reserialized DB JSON.
        for i,(r,s) in enumerate(zip(all_reads,sources)):
            ordinals[s['input_id']] += 1
            valid = valid and (r['read_index'] == i and s['query_name'] == f'q{i}'
                     and s['record_ordinal'] == ordinals[s['input_id']]
                     and s['sequence_sha256'] == evidence['query_sha256'][i])
    except (KeyError, TypeError, IndexError, ValueError):
        valid = False
    if not valid:
        raise HTTPException(409, 'alignment_source_mapping_invalid')


@router.post('', status_code=201)
def submit(revision_id: str, request: Submission):
    with closing(connect()) as db:
        _revision(db, revision_id)
    try:
        job = enqueue_alignment(revision_id, request.input_ids, request.data_type,
                                request.idempotency_key, request.max_attempts)
    except QueueConflict as exc:
        raise HTTPException(409, str(exc)) from exc
    except (InputValidationError, ValueError) as exc:
        raise HTTPException(422, str(exc)) from exc
    with closing(connect()) as db:
        return _summary(_job(db, revision_id, job))


@router.get('')
def listing(revision_id: str, limit: int=Query(5,ge=1,le=100), offset: int=Query(0,ge=0)):
    with closing(connect()) as db:
        db.execute('BEGIN');_revision(db, revision_id)
        query = 'FROM analysis_jobs j JOIN queued_jobs q ON q.job_id=j.id WHERE j.sequence_revision_id=? AND q.adapter=?'
        args = (revision_id, ADAPTER)
        total = db.execute('SELECT count(*) '+query, args).fetchone()[0]
        ids = db.execute('SELECT j.id '+query+' ORDER BY j.created_at DESC,j.id DESC LIMIT ? OFFSET ?',
                         (*args, limit, offset)).fetchall()
        return _page([_summary(_job(db,revision_id,r['id'])) for r in ids], total, limit, offset)


@router.get('/{job_id}')
def detail(revision_id: str, job_id: str, limit: int=Query(5,ge=1,le=100), offset: int=Query(0,ge=0)):
    with closing(connect()) as db:
        db.execute('BEGIN');row = _job(db, revision_id, job_id)
        response = _summary(row)
        manifest = json.loads(row['input_manifest_json'])
        response['reference'] = manifest['reference']
        response['input_identities'] = manifest['inputs']['fastq_inputs']
        response['tool'] = manifest['inputs']['tool']
        response['input_labels'] = {}
        for item in response['input_identities']:
            label = db.execute('SELECT original_filename FROM fastq_inputs WHERE id=? AND sequence_revision_id=?',
                               (item['id'],revision_id)).fetchone()
            response['input_labels'][item['id']] = label[0] if label else None
        response['result'] = None
        if row['status'] == 'succeeded':
            evidence = _published(row)
            response['result'] = {k:evidence[k] for k in RESULT_FIELDS}
            response['result']['read_count'] = len(evidence['reads'])
        total = db.execute('SELECT count(*) FROM job_attempts WHERE job_id=?',(job_id,)).fetchone()[0]
        attempts = [dict(r) for r in db.execute('''SELECT number,status,started_at,completed_at,error_detail
                    FROM job_attempts WHERE job_id=? ORDER BY number DESC LIMIT ? OFFSET ?''',
                    (job_id,limit,offset))]
        response['attempts'] = _page(attempts,total,limit,offset)
        return response


@router.get('/{job_id}/reads')
def reads(revision_id: str, job_id: str, limit: int=Query(20,ge=1,le=100), offset: int=Query(0,ge=0),
          start: Optional[int]=Query(None,ge=0), end: Optional[int]=Query(None,ge=0),
          relation: Literal['paired','deletion','either']='either'):
    with closing(connect()) as db:
        db.execute('BEGIN');evidence = _published(_job(db,revision_id,job_id))
    _validate_sources(evidence)
    all_reads, sources = evidence['reads'], evidence['sources']
    selected = None
    if (start is None) != (end is None) or (start is None and relation != 'either'):
        raise HTTPException(422,'both_reference_boundaries_required')
    if start is not None:
        from .alignment_coverage import select_region, CoverageError
        try:
            selected = select_region(evidence,start,end,relation)
        except CoverageError as exc:
            raise HTTPException(422 if str(exc)=='invalid_reference_region' else 409,str(exc)) from exc
    indices = list(range(len(all_reads))) if selected is None else list(selected)
    items = [{**{k:all_reads[i][k] for k in READ_FIELDS},
              'source':{k:sources[i][k] for k in ('query_name','input_id','record_ordinal','sequence_sha256')},
              **({'region_evidence':selected[i]} if selected is not None else {})}
             for i in indices[offset:offset+limit]]
    return {**_page(items,len(indices),limit,offset), 'scope':SCOPE,
            'unfiltered_total':len(all_reads),
            'region':None if selected is None else {'start':start,'end':end,'relation':relation},
            'coordinate_convention':evidence['coordinate_convention'],
            'base_quality_used':evidence['base_quality_used'],
            'whole_reference_verified':evidence['whole_reference_verified']}


@router.post('/{job_id}/cancel')
def cancellation(revision_id: str, job_id: str):
    with closing(connect()) as db:
        _job(db,revision_id,job_id)
    try:
        cancel(job_id)
    except QueueConflict as exc:
        raise HTTPException(409,str(exc)) from exc
    with closing(connect()) as db:
        return _summary(_job(db,revision_id,job_id))


@router.get('/{job_id}/coverage')
def coverage(revision_id: str, job_id: str,
             kind: Literal['unpaired','ambiguous_only','deletion']='unpaired',
             limit: int=Query(20,ge=1,le=100), offset: int=Query(0,ge=0)):
    from .alignment_coverage import summarize, regions, CoverageError
    with closing(connect()) as db:
        db.execute('BEGIN');evidence = _published(_job(db,revision_id,job_id))
    _validate_sources(evidence)
    try:
        summary = summarize(evidence)
        intervals = regions(summary,kind)
    except CoverageError as exc:
        raise HTTPException(409,str(exc)) from exc
    summary.pop('segments')
    return {'job_id':job_id,'attempt_number':evidence['attempt_number'],
            'reference_sha256':evidence['reference_sha256'],'summary':summary,
            'regions':{**_page([{'start':a,'end':b,'length_bp':b-a} for a,b in intervals[offset:offset+limit]],
                              len(intervals),limit,offset),'kind':kind}}
