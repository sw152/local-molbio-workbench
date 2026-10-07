"""Revision-scoped input-check operations. Worker lease credentials stay internal."""
from contextlib import closing
import json
from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, ConfigDict, Field
from .database import connect
from .input_validation import ADAPTER, InputValidationError, enqueue_check
from .job_queue import QueueConflict, cancel

router = APIRouter(prefix='/api/revisions/{revision_id}/input-checks')


class Submission(BaseModel):
    model_config = ConfigDict(extra='forbid')
    read_ids: list[str] = Field(min_length=1, max_length=100)
    idempotency_key: str = Field(min_length=1, max_length=200)
    max_attempts: int = Field(default=3, ge=1, le=10, strict=True)


def _revision(db, revision):
    if db.execute('SELECT 1 FROM sequence_revisions WHERE id=?',(revision,)).fetchone() is None:
        raise HTTPException(404,'Sequence revision not found')


def _job(db, revision, job):
    row=db.execute('''SELECT j.id,j.status,j.created_at,j.started_at,j.completed_at,j.error_detail,
                     j.input_manifest_json,j.result_summary_json,q.attempt_count,q.max_attempts,q.lease_expires
                     FROM analysis_jobs j JOIN queued_jobs q ON q.job_id=j.id
                     WHERE j.id=? AND j.sequence_revision_id=? AND q.adapter=?''',(job,revision,ADAPTER)).fetchone()
    if row is None:raise HTTPException(404,'Input check not found for this revision')
    return row


def _summary(row):
    result=dict(row)
    manifest=json.loads(result.pop('input_manifest_json'))
    result.pop('result_summary_json')
    result['read_count']=len(manifest['inputs']['reads'])
    result['scope']='registered_input_integrity_only'
    return result


@router.post('',status_code=201)
def submit(revision_id: str, request: Submission):
    with closing(connect()) as db:_revision(db,revision_id)
    try:job=enqueue_check(revision_id,request.read_ids,request.idempotency_key,request.max_attempts)
    except QueueConflict as exc:raise HTTPException(409,str(exc)) from exc
    except InputValidationError as exc:raise HTTPException(422,str(exc)) from exc
    except ValueError as exc:raise HTTPException(422,str(exc)) from exc
    with closing(connect()) as db:return _summary(_job(db,revision_id,job))


@router.get('')
def listing(revision_id: str,limit: int=Query(5,ge=1,le=100),offset: int=Query(0,ge=0)):
    with closing(connect()) as db:
        db.execute('BEGIN');_revision(db,revision_id)
        query='FROM analysis_jobs j JOIN queued_jobs q ON q.job_id=j.id WHERE j.sequence_revision_id=? AND q.adapter=?'
        args=(revision_id,ADAPTER)
        total=db.execute('SELECT count(*) '+query,args).fetchone()[0]
        ids=db.execute('SELECT j.id '+query+' ORDER BY j.created_at DESC,j.id DESC LIMIT ? OFFSET ?',(*args,limit,offset)).fetchall()
        items=[_summary(_job(db,revision_id,row['id'])) for row in ids]
    return {'items':items,'total':total,'offset':offset,'limit':limit,'has_more':offset+len(items)<total}


@router.get('/{job_id}')
def detail(revision_id: str,job_id: str,limit: int=Query(5,ge=1,le=100),offset: int=Query(0,ge=0)):
    with closing(connect()) as db:
        db.execute('BEGIN');row=_job(db,revision_id,job_id)
        result=_summary(row)
        result['inputs']=json.loads(row['input_manifest_json'])
        result['result']=json.loads(row['result_summary_json'])
        result['read_labels']={}
        for item in result['inputs']['inputs']['reads']:
            read=db.execute('SELECT original_filename FROM sequencing_reads WHERE id=? AND sequence_revision_id=?',(item['id'],revision_id)).fetchone()
            result['read_labels'][item['id']]=read['original_filename'] if read else None
        total=db.execute('SELECT count(*) FROM job_attempts WHERE job_id=?',(job_id,)).fetchone()[0]
        attempts=[dict(r) for r in db.execute('''SELECT number,worker_id,status,started_at,completed_at,error_detail
                     FROM job_attempts WHERE job_id=? ORDER BY number DESC LIMIT ? OFFSET ?''',(job_id,limit,offset))]
        result['attempts']={'items':attempts,'total':total,'offset':offset,'limit':limit,'has_more':offset+len(attempts)<total}
    return result


@router.post('/{job_id}/cancel')
def cancellation(revision_id: str,job_id: str):
    with closing(connect()) as db:_job(db,revision_id,job_id)
    try:cancel(job_id)
    except QueueConflict as exc:raise HTTPException(409,str(exc)) from exc
    with closing(connect()) as db:return _summary(_job(db,revision_id,job_id))
