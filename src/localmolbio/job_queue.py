"""Durable single-host queue primitives. No analysis adapter or daemon is started here.

Leases provide at-least-once execution with fenced result publication, not exactly-once
external side effects. Workers must stage artifacts separately for each attempt.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
from hashlib import sha256
import json
import time
from uuid import uuid4

from .database import connect

SCHEMA = """
CREATE TABLE IF NOT EXISTS queued_jobs (
    job_id TEXT PRIMARY KEY REFERENCES analysis_jobs(id) ON DELETE RESTRICT,
    idempotency_key TEXT NOT NULL UNIQUE,
    request_sha256 TEXT NOT NULL,
    adapter TEXT NOT NULL,
    max_attempts INTEGER NOT NULL CHECK(max_attempts BETWEEN 1 AND 10),
    attempt_count INTEGER NOT NULL DEFAULT 0,
    lease_token TEXT,
    lease_expires REAL
);
CREATE TABLE IF NOT EXISTS job_attempts (
    id TEXT PRIMARY KEY,
    job_id TEXT NOT NULL REFERENCES queued_jobs(job_id) ON DELETE RESTRICT,
    number INTEGER NOT NULL,
    worker_id TEXT NOT NULL,
    status TEXT NOT NULL CHECK(status IN ('running','succeeded','failed','expired','cancelled')),
    started_at TEXT NOT NULL,
    completed_at TEXT,
    error_detail TEXT,
    UNIQUE(job_id,number)
);
"""

class QueueConflict(ValueError):
    pass


def _json(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False)


def _stamp(now):
    return datetime.fromtimestamp(now,timezone.utc).isoformat()


def _name(value):
    if not isinstance(value,str) or not value.strip() or len(value)>200:
        raise ValueError('Expected a nonempty identifier of at most 200 characters')
    return value


@contextmanager
def _transaction():
    db=connect()
    try:
        db.execute('BEGIN IMMEDIATE')
        yield db
        db.commit()
    except BaseException:
        db.rollback()
        raise
    finally:
        db.close()


def _audit(db,job,action,payload,now):
    db.execute('INSERT INTO audit_events VALUES (?,?,?,?,?,?)',
               (str(uuid4()),'analysis-job',job,action,_json(payload),_stamp(now)))


def enqueue(revision_id, idempotency_key, adapter, input_manifest, parameters=None, max_attempts=3):
    """Store an immutable request; a reused key must describe exactly the same request."""
    _name(idempotency_key);_name(adapter)
    if type(max_attempts) is not int or not 1<=max_attempts<=10:
        raise ValueError('max_attempts must be an integer from 1 to 10')
    if not isinstance(input_manifest,dict) or (parameters is not None and not isinstance(parameters,dict)):
        raise ValueError('Manifest and parameters must be objects')
    parameters={} if parameters is None else parameters
    with _transaction() as db:
        reference=db.execute('SELECT id,sequence_sha256,topology FROM sequence_revisions WHERE id=?',(revision_id,)).fetchone()
        if reference is None:raise ValueError('Sequence revision not found')
        manifest={'reference':dict(reference),'inputs':input_manifest,'adapter':adapter}
        request=_json({'manifest':manifest,'parameters':parameters,'max_attempts':max_attempts})
        digest=sha256(request.encode()).hexdigest()
        existing=db.execute('SELECT job_id,request_sha256 FROM queued_jobs WHERE idempotency_key=?',(idempotency_key,)).fetchone()
        if existing:
            if existing['request_sha256']!=digest:raise QueueConflict('Idempotency key already belongs to a different request')
            return existing['job_id']
        now=time.time();job=str(uuid4())
        db.execute('''INSERT INTO analysis_jobs (id,sequence_revision_id,job_kind,status,parameters_json,input_manifest_json,created_at)
                      VALUES (?,?,'plasmid-verification','queued',?,?,?)''',
                   (job,revision_id,_json(parameters),_json(manifest),_stamp(now)))
        db.execute('INSERT INTO queued_jobs (job_id,idempotency_key,request_sha256,adapter,max_attempts) VALUES (?,?,?,?,?)',
                   (job,idempotency_key,digest,adapter,max_attempts))
        _audit(db,job,'queued',{'adapter':adapter,'request_sha256':digest},now)
        return job


def _recover(db,now):
    rows=db.execute('''SELECT q.* FROM queued_jobs q JOIN analysis_jobs j ON j.id=q.job_id
                       WHERE j.status='running' AND q.lease_expires<=?''',(now,)).fetchall()
    for row in rows:
        retry=row['attempt_count']<row['max_attempts'];job=row['job_id']
        db.execute("UPDATE job_attempts SET status='expired',completed_at=?,error_detail='Worker lease expired' WHERE id=?",
                   (_stamp(now),row['lease_token']))
        db.execute('UPDATE analysis_jobs SET status=?,completed_at=?,error_detail=? WHERE id=?',
                   ('queued' if retry else 'failed',None if retry else _stamp(now),'Worker lease expired',job))
        db.execute('UPDATE queued_jobs SET lease_token=NULL,lease_expires=NULL WHERE job_id=?',(job,))
        _audit(db,job,'lease-expired',{'attempt':row['attempt_count'],'requeued':retry},now)
    return len(rows)


def recover_expired():
    with _transaction() as db:return _recover(db,time.time())


def _duration(seconds):
    if type(seconds) is not int or not 1<=seconds<=3600:
        raise ValueError('Lease duration must be 1–3600 integer seconds')


def claim(worker_id, adapter, lease_seconds=60):
    """Claim one compatible job, recovering expired leases inside the same lock."""
    _name(worker_id);_name(adapter);_duration(lease_seconds)
    with _transaction() as db:
        now=time.time();_recover(db,now)
        row=db.execute('''SELECT j.*,q.attempt_count FROM analysis_jobs j JOIN queued_jobs q ON q.job_id=j.id
                          WHERE j.status='queued' AND q.adapter=? AND q.attempt_count<q.max_attempts
                          ORDER BY j.created_at,j.id LIMIT 1''',(adapter,)).fetchone()
        if row is None:return None
        job=row['id'];number=row['attempt_count']+1;token=str(uuid4());expires=now+lease_seconds
        db.execute('UPDATE queued_jobs SET attempt_count=?,lease_token=?,lease_expires=? WHERE job_id=?',(number,token,expires,job))
        db.execute("UPDATE analysis_jobs SET status='running',started_at=COALESCE(started_at,?),completed_at=NULL,error_detail=NULL WHERE id=?",(_stamp(now),job))
        db.execute("INSERT INTO job_attempts (id,job_id,number,worker_id,status,started_at) VALUES (?,?,?,?,'running',?)",(token,job,number,worker_id,_stamp(now)))
        _audit(db,job,'claimed',{'attempt':number,'worker_id':worker_id},now)
        return {'job_id':job,'attempt':number,'lease_token':token,'lease_expires':expires,
                'adapter':adapter,'input_manifest':json.loads(row['input_manifest_json']),'parameters':json.loads(row['parameters_json'])}


def _owned(db,job,token,now):
    row=db.execute('''SELECT q.*,j.status FROM queued_jobs q JOIN analysis_jobs j ON j.id=q.job_id WHERE q.job_id=?''',(job,)).fetchone()
    if row is None or row['status']!='running' or row['lease_token']!=token or row['lease_expires']<=now:
        raise QueueConflict('Lease is missing, expired or superseded')
    return row


def renew(job_id, lease_token, lease_seconds=60):
    _duration(lease_seconds)
    with _transaction() as db:
        now=time.time();_owned(db,job_id,lease_token,now)
        expires=now+lease_seconds
        db.execute('UPDATE queued_jobs SET lease_expires=? WHERE job_id=?',(expires,job_id))
        return expires


def finish(job_id, lease_token, result):
    if not isinstance(result,dict):raise ValueError('Result must be an object')
    encoded=_json(result)
    with _transaction() as db:
        now=time.time();_owned(db,job_id,lease_token,now)
        db.execute("UPDATE analysis_jobs SET status='succeeded',completed_at=?,result_summary_json=?,error_detail=NULL WHERE id=?",(_stamp(now),encoded,job_id))
        db.execute("UPDATE job_attempts SET status='succeeded',completed_at=? WHERE id=?",(_stamp(now),lease_token))
        db.execute('UPDATE queued_jobs SET lease_token=NULL,lease_expires=NULL WHERE job_id=?',(job_id,))
        _audit(db,job_id,'succeeded',{'attempt_id':lease_token},now)


def fail(job_id, lease_token, error, retryable=False):
    if not isinstance(error,str) or not error or len(error)>4000:raise ValueError('Error must be 1–4000 characters')
    if type(retryable) is not bool:raise ValueError('retryable must be boolean')
    with _transaction() as db:
        now=time.time();row=_owned(db,job_id,lease_token,now)
        retry=retryable and row['attempt_count']<row['max_attempts']
        db.execute('UPDATE analysis_jobs SET status=?,completed_at=?,error_detail=? WHERE id=?',('queued' if retry else 'failed',None if retry else _stamp(now),error,job_id))
        db.execute("UPDATE job_attempts SET status='failed',completed_at=?,error_detail=? WHERE id=?",(_stamp(now),error,lease_token))
        db.execute('UPDATE queued_jobs SET lease_token=NULL,lease_expires=NULL WHERE job_id=?',(job_id,))
        _audit(db,job_id,'attempt-failed',{'attempt_id':lease_token,'requeued':retry},now)


def cancel(job_id):
    with _transaction() as db:
        row=db.execute('SELECT q.lease_token,j.status FROM queued_jobs q JOIN analysis_jobs j ON j.id=q.job_id WHERE q.job_id=?',(job_id,)).fetchone()
        if row is None:raise ValueError('Queued job not found')
        if row['status']=='cancelled':return
        if row['status'] not in {'queued','running'}:raise QueueConflict('Completed jobs cannot be cancelled')
        now=time.time()
        db.execute("UPDATE analysis_jobs SET status='cancelled',completed_at=? WHERE id=?",(_stamp(now),job_id))
        if row['lease_token']:
            db.execute("UPDATE job_attempts SET status='cancelled',completed_at=? WHERE id=?",(_stamp(now),row['lease_token']))
        db.execute('UPDATE queued_jobs SET lease_token=NULL,lease_expires=NULL WHERE job_id=?',(job_id,))
        _audit(db,job_id,'cancelled',{},now)
