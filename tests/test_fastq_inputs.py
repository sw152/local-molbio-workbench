"""FASTQ syntax, quality boundaries, compressed integrity and atomic registration."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import gzip
import io
import json
from pathlib import Path
import random

import pytest
from Bio import SeqIO
from localmolbio.fastq import FastqError, FastqLimits, inspect_fastq
from localmolbio.fastq_inputs import FastqConflict, register_fastq
from localmolbio.database import connect, initialise
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision

DATA=b'@first\nACGN\n+first\n!45?\n@second\nry\n+\n@+\n'


def inspect(tmp_path,data=DATA,compressed=False,**kwargs):
    path=tmp_path/'input';path.write_bytes(gzip.compress(data,mtime=0) if compressed else data)
    return inspect_fastq(path,compressed,'phred33',**kwargs)


@pytest.mark.parametrize('compressed',[False,True])
def test_quality_and_ambiguity_statistics_are_exact(tmp_path,compressed):
    s=inspect(tmp_path,compressed=compressed)
    assert s['records']==2 and s['bases']==6
    assert (s['min_read_length'],s['max_read_length'])==(2,4)
    assert s['n_bases']==1 and s['ambiguous_bases']==3
    assert (s['min_phred'],s['max_phred'],s['q20_bases'],s['q30_bases'])==(0,31,3,2)
    assert s['mean_phred']==pytest.approx((0+19+20+30+31+10)/6)
    assert not s['analysis_performed'] and not s['whole_reference_verified']


def test_wrapped_crlf_quality_at_and_plus_match_biopython(tmp_path):
    rng=random.Random(902)
    chunks=[]
    for i in range(20):
        sequence=''.join(rng.choice('ACGTN') for _ in range(30))
        quality='@+'+''.join(chr(rng.randrange(33,127)) for _ in range(28))
        chunks.append(f'@r{i}\r\n{sequence[:15]}\r\n{sequence[15:]}\r\n+r{i}\r\n{quality[:1]}\r\n{quality[1:]}\r\n')
    text=''.join(chunks).rstrip('\r\n')
    records=list(SeqIO.parse(io.StringIO(text),'fastq'))
    s=inspect(tmp_path,text.encode())
    expected=[q for r in records for q in r.letter_annotations['phred_quality']]
    assert s['records']==len(records)==20 and s['bases']==600
    assert s['mean_phred']==pytest.approx(sum(expected)/len(expected))
    assert s['q20_bases']==sum(q>=20 for q in expected)


@pytest.mark.parametrize('data,reason',[
    (b'','empty_fastq'),(b'\n','invalid_fastq_header'),(b'@\nAC\n+\nII\n','invalid_fastq_header'),
    (b'@x\nAC\n','truncated_sequence'),(b'@x\nAC\n+y\nII\n','repeated_header_mismatch'),
    (b'@x\n+\n','empty_read'),(b'@x\nAU\n+\nII\n','invalid_dna_sequence'),
    (b'@x\nAC G\n+\nIIII\n','invalid_dna_sequence'),(b'@x\nAC\n+\nI','truncated_quality'),
    (b'@x\nAC\n+\nIII','sequence_quality_length_mismatch'),
    (b'@x\nAC\n+\n I','invalid_phred33_character'),(b'@x\nAC\n+\nI\x7f','invalid_phred33_character'),
    (DATA+b'trailing','invalid_fastq_header'),
])
def test_invalid_records_are_rejected_without_partial_results(tmp_path,data,reason):
    with pytest.raises(FastqError,match=reason):inspect(tmp_path,data)


def test_gzip_crc_truncation_concatenation_and_filename_mismatch(tmp_path):
    path=tmp_path/'input'
    valid=gzip.compress(DATA,mtime=0)
    for broken in [valid[:-4],valid[:-8]+b'\0'*8,valid+b'not-gzip']:
        path.write_bytes(broken)
        with pytest.raises(FastqError,match='invalid_or_truncated_gzip'):inspect_fastq(path,True,'phred33')
    path.write_bytes(valid+valid)
    assert inspect_fastq(path,True,'phred33')['records']==4
    with pytest.raises(FastqError,match='compression_does_not_match_filename'):inspect_fastq(path,False,'phred33')
    path.write_bytes(DATA)
    with pytest.raises(FastqError,match='compression_does_not_match_filename'):inspect_fastq(path,True,'phred33')
    with pytest.raises(FastqError,match='explicit_phred33'):inspect_fastq(path,False,'auto')


@pytest.mark.parametrize('limits,reason',[
    (FastqLimits(raw_bytes=5),'raw_size_limit_exceeded'),
    (FastqLimits(decoded_bytes=10),'decoded_size_limit_exceeded'),
    (FastqLimits(records=1),'record_count_limit_exceeded'),
    (FastqLimits(record_bases=2),'line_size_limit_exceeded'),
])
def test_limits_reject_without_silent_truncation(tmp_path,limits,reason):
    with pytest.raises(FastqError,match=reason):inspect(tmp_path,compressed=True,limits=limits)


def test_multiline_record_limit_and_low_quality_are_not_trimmed(tmp_path):
    with pytest.raises(FastqError,match='record_size_limit_exceeded'):
        inspect(tmp_path,b'@x\nAC\nGT\n+\n!!!!\n',limits=FastqLimits(record_bases=3))
    s=inspect(tmp_path,b'@x\nACGT\n+\n!!!!\n')
    assert s['bases']==4 and s['q20_bases']==0
    with pytest.raises(ValueError):FastqLimits(records=True)


def endpoint(revision):return f'/api/revisions/{revision}/fastq-inputs'
def upload(client,revision,data=DATA,name='reads.fastq',encoding='phred33'):
    return client.post(endpoint(revision),data={'quality_encoding':encoding},files={'file':(name,data)})


def test_api_preserves_original_bytes_and_keeps_fastq_out_of_sanger(ready):
    client,_,revision=ready;blob=gzip.compress(DATA,mtime=0)
    result=upload(client,revision,blob,'reads.fastq.gz')
    assert result.status_code==201,result.text
    value=result.json();assert value['file_sha256']==sha256(blob).hexdigest()
    with connect() as db:
        path=Path(db.execute('SELECT storage_path FROM fastq_inputs WHERE id=?',(value['id'],)).fetchone()[0])
        assert db.execute("SELECT count(*) FROM audit_events WHERE object_type='fastq-input'").fetchone()[0]==1
    assert path.read_bytes()==blob and not list((path.parent/'.incoming').iterdir())
    assert 'storage_path' not in result.text and 'ACGN' not in result.text
    assert len(client.get(f'/api/sequence-revisions/{revision}/sanger-reads').json())==1
    assert client.get(endpoint(revision)).json()['items']==[value]
    assert upload(client,revision,blob,'same.fq.gz').status_code==409
    other=sibling_revision(revision)
    assert client.get(endpoint(other)).json()['total']==0
    assert upload(client,other,blob,'other.fq.gz').status_code==201
    with connect() as db:assert db.execute('SELECT count(DISTINCT storage_path) FROM fastq_inputs').fetchone()[0]==1
    initialise();assert client.get(endpoint(revision)).json()['total']==1


def test_rejected_uploads_leave_no_records_or_staging_files(ready):
    client,_,revision=ready
    for data,name,encoding in [(b'bad','bad.fastq','phred33'),(DATA,'bad.zip','phred33'),(DATA,'x.fastq','auto'),(gzip.compress(DATA),'wrong.fastq','phred33')]:
        assert upload(client,revision,data,name,encoding).status_code==422
    assert client.post(endpoint(revision),files={'file':('reads.fq',DATA)}).status_code==422
    assert upload(client,'missing').status_code==404
    with connect() as db:assert db.execute('SELECT count(*) FROM fastq_inputs').fetchone()[0]==0
    from localmolbio.config import data_dir
    directory=data_dir()/'fastq'
    assert not list(directory.glob('*.fastq*')) and not list((directory/'.incoming').iterdir())
    with pytest.raises(FastqError,match='raw_size_limit_exceeded'):
        register_fastq(revision,'large.fq',io.BytesIO(DATA),'phred33',FastqLimits(raw_bytes=2))
    assert not list((directory/'.incoming').iterdir())


def test_registration_is_concurrent_safe_and_existing_bytes_never_overwritten(ready):
    _,_,revision=ready
    def register(_):
        try:return register_fastq(revision,'same.fq',io.BytesIO(DATA),'phred33')['id']
        except FastqConflict:return None
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(register,range(4)))
    assert sum(r is not None for r in results)==1
    with connect() as db:path=Path(db.execute('SELECT storage_path FROM fastq_inputs').fetchone()[0])
    path.write_bytes(b'externally modified')
    other=sibling_revision(revision)
    with pytest.raises(FastqConflict,match='existing_fastq_storage_integrity_mismatch'):
        register_fastq(other,'same.fq',io.BytesIO(DATA),'phred33')
    assert path.read_bytes()==b'externally modified'
    with connect() as db:assert db.execute('SELECT count(*) FROM fastq_inputs').fetchone()[0]==1


def test_pagination_returns_metadata_only(ready):
    client,_,revision=ready
    for i in range(3):assert upload(client,revision,DATA.replace(b'first',f'read{i}'.encode()),f'{i}.fq').status_code==201
    first=client.get(endpoint(revision)+'?limit=2').json();second=client.get(endpoint(revision)+'?limit=2&offset=2').json()
    assert first['has_more'] and not second['has_more'] and first['total']==3
    assert len({r['id'] for r in first['items']+second['items']})==3
    assert client.get(endpoint(revision)+'?limit=0').status_code==422
    assert client.get(endpoint('missing')).status_code==404
