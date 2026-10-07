from io import BytesIO, StringIO
from zipfile import ZipFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastapi.testclient import TestClient
from localmolbio.api import app
from localmolbio.database import connect


def test_library_pagination_reaches_every_record_and_searches_literal_names(tmp_path,monkeypatch):
    monkeypatch.setenv('MOLBIO_DATA_DIR',str(tmp_path/'runtime'))
    archive=BytesIO()
    with ZipFile(archive,'w') as z:
        for i in range(205):
            name=f'record-{i:03d}'
            record=SeqRecord(Seq('ACGT'*5),id=name,name=name,annotations={'molecule_type':'DNA'})
            stream=StringIO();SeqIO.write(record,stream,'genbank');z.writestr(name+'.gb',stream.getvalue())
    with TestClient(app) as client:
        assert client.post('/api/imports/benchling',files={'file':('library.zip',archive.getvalue())}).status_code==201
        all_ids=[]
        for offset in range(0,205,24):
            response=client.get('/api/sequence-library',params={'offset':offset,'limit':24})
            assert response.status_code==200
            page=response.json()
            assert page['total']==205
            assert page['offset']==offset
            assert page['has_more']==(offset+24<205)
            all_ids.extend(item['id'] for item in page['items'])
        assert len(all_ids)==len(set(all_ids))==205
        assert client.get('/api/sequence-library?offset=205').json()['items']==[]
        assert len(client.get('/api/sequences').json())==200  # old array interface remains compatible
        with connect() as db:
            # Equal names must still have stable page boundaries.
            db.execute("UPDATE sequences SET display_name='same name'")
        first=client.get('/api/sequence-library?limit=24').json()['items']
        assert [r['id'] for r in first]==all_ids[:24]
        with connect() as db:
            db.execute('UPDATE sequences SET display_name=? WHERE id=?',('100%_exact!',all_ids[0]))
        for query in ['%','_','!','100%_exact!']:
            page=client.get('/api/sequence-library',params={'query':query}).json()
            assert page['total']==1 and page['items'][0]['id']==all_ids[0]
        assert client.get('/api/sequence-library',params={'query':'record-204.gb'}).json()['total']==1
        assert client.get('/api/sequence-library',params={'query':'not present'}).json()['total']==0
        for params in [{'offset':-1},{'limit':0},{'limit':101},{'query':'a'*501}]:
            assert client.get('/api/sequence-library',params=params).status_code==422
