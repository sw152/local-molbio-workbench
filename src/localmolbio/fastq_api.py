from contextlib import closing
import json
from typing import Literal
from fastapi import APIRouter, File, Form, HTTPException, Query, UploadFile
from .database import connect
from .fastq import FastqError
from .fastq_inputs import FastqConflict, register_fastq

router=APIRouter(prefix='/api/revisions/{revision_id}/fastq-inputs')


@router.post('',status_code=201)
def upload(revision_id: str, quality_encoding: Literal['phred33']=Form(...), file: UploadFile=File(...)):
    try:return register_fastq(revision_id,file.filename,file.file,quality_encoding)
    except FastqConflict as exc:raise HTTPException(409,str(exc)) from exc
    except FastqError as exc:raise HTTPException(404 if str(exc)=='reference_missing' else 422,str(exc)) from exc


@router.get('')
def listing(revision_id: str,limit: int=Query(10,ge=1,le=100),offset: int=Query(0,ge=0)):
    with closing(connect()) as db:
        db.execute('BEGIN')
        if db.execute('SELECT 1 FROM sequence_revisions WHERE id=?',(revision_id,)).fetchone() is None:
            raise HTTPException(404,'Sequence revision not found')
        total=db.execute('SELECT count(*) FROM fastq_inputs WHERE sequence_revision_id=?',(revision_id,)).fetchone()[0]
        rows=db.execute('''SELECT id,sequence_revision_id,original_filename,file_sha256,size_bytes,compression,quality_encoding,summary_json,created_at
                           FROM fastq_inputs WHERE sequence_revision_id=? ORDER BY created_at DESC,id DESC LIMIT ? OFFSET ?''',(revision_id,limit,offset)).fetchall()
        items=[]
        for row in rows:
            item=dict(row);item['summary']=json.loads(item.pop('summary_json'));items.append(item)
    return {'items':items,'total':total,'offset':offset,'limit':limit,'has_more':offset+len(items)<total}
