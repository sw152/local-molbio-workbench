from __future__ import annotations

from io import BytesIO
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi.testclient import TestClient

from localmolbio.api import app
from localmolbio.database import connect
from localmolbio.sanger import SangerRead


def test_import_endpoint_exposes_annotation_review_status(tmp_path, monkeypatch) -> None:
    archive = BytesIO()
    record = """LOCUS       wrapped                 8 bp    ds-DNA     circular     01-JAN-2026
DEFINITION  synthetic test record.
FEATURES             Location/Qualifiers
     source          1..8
                     /organism=\"synthetic construct\"
     misc_feature    6..3
ORIGIN
        1 atgcatgc
//
"""
    with ZipFile(archive, "w", ZIP_DEFLATED) as zipper:
        zipper.writestr("wrapped.gb", record)

    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    with TestClient(app) as client:
        response = client.post(
            "/api/imports/benchling",
            files={"file": ("benchling.zip", archive.getvalue(), "application/zip")},
        )
        assert response.status_code == 201
        assert response.json()["parse_warning_records"] == 1

        sequences = client.get("/api/sequences")
        assert sequences.status_code == 200
        assert sequences.json()[0]["parse_warning_count"] == 1

        revisions = client.get(f"/api/sequences/{sequences.json()[0]['id']}/revisions")
        assert revisions.status_code == 200
        assert revisions.json()[0]["revision_number"] == 1
        assert revisions.json()[0]["source_kind"] == "import"

        detail = client.get(f"/api/sequences/{sequences.json()[0]['id']}")
        assert detail.status_code == 200
        assert "spans the origin" in detail.json()["parse_warnings_json"]

        with connect() as database:
            database.execute("UPDATE sequences SET features_json = '[]'")
        map_response = client.get(f"/api/sequences/{sequences.json()[0]['id']}/map")
        assert map_response.status_code == 200
        assert map_response.json()["requires_annotation_review"] is True
        assert map_response.json()["features"][0]["renderable"] is True


def test_sanger_upload_persists_parsed_read(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    monkeypatch.setattr(
        "localmolbio.api.parse_ab1",
        lambda _: SangerRead(
            sequence="ATGCATGC",
            quality_summary={"mean_phred": 35.0, "q20_fraction": 1.0},
            parser_metadata={"SMPL1": "synthetic"},
        ),
    )
    with TestClient(app) as client:
        revision_id = client.get("/api/sequences").json()
        assert revision_id == []
        # Seed a minimal sequence through the existing import route.
        archive = BytesIO()
        with ZipFile(archive, "w", ZIP_DEFLATED) as zipper:
            zipper.writestr(
                "sequence.gb",
                """LOCUS       one                      8 bp    ds-DNA     circular     01-JAN-2026
DEFINITION  synthetic test record.
FEATURES             Location/Qualifiers
     source          1..8
ORIGIN
        1 atgcatgc
//
""",
            )
        imported = client.post(
            "/api/imports/benchling",
            files={"file": ("benchling.zip", archive.getvalue(), "application/zip")},
        )
        sequence_id = client.get("/api/sequences").json()[0]["id"]
        revision_id = client.get(f"/api/sequences/{sequence_id}/revisions").json()[0]["id"]
        uploaded = client.post(
            "/api/sanger-reads",
            data={"sequence_revision_id": revision_id, "direction": "forward"},
            files={"file": ("trace.ab1", b"synthetic-ab1", "application/octet-stream")},
        )
    assert imported.status_code == 201
    assert uploaded.status_code == 201
    assert uploaded.json()["quality_summary"]["mean_phred"] == 35.0


def test_sanger_evidence_survives_reload_and_records_inputs(tmp_path, monkeypatch):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from io import StringIO
    import json

    from abif_fixture import synthetic_ab1
    reference = "ATCGGATCAGGTACGTTAGCTACGTGGTACCACTGATCCGAGTACGTAGCTACGACATCG"
    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    monkeypatch.setattr("localmolbio.api.parse_ab1", lambda _: SangerRead(
        reference[-30:], {"mean_phred": 40, "q20_fraction": 1}, {}, [40] * 30,
    ))
    record = SeqRecord(Seq(reference), id="synthetic", name="synthetic", annotations={"molecule_type": "DNA", "topology": "linear"})
    handle = StringIO()
    SeqIO.write(record, handle, "genbank")
    archive = BytesIO()
    with ZipFile(archive, "w") as zipper:
        zipper.writestr("synthetic.gb", handle.getvalue())
    with TestClient(app) as client:
        assert client.post("/api/imports/benchling", files={"file": ("synthetic.zip", archive.getvalue())}).status_code == 201
        sequence_id = client.get("/api/sequences").json()[0]["id"]
        revision = client.get(f"/api/sequences/{sequence_id}/revisions").json()[0]["id"]
        uploaded = client.post("/api/sanger-reads", data={"sequence_revision_id": revision}, files={"file": ("synthetic.ab1", synthetic_ab1(reference[-30:], [40] * 30, trace=True))})
        assert uploaded.status_code == 201
        result = client.post("/api/sanger-verifications", json={"sequencing_read_id": uploaded.json()["id"]})
        assert result.status_code == 201, result.text
        evidence = result.json()["evidence"]
        assert evidence["direction"] == "forward"
        assert evidence["review_flags"] == []
        assert result.json()["reference_end"] == len(reference)
        assert evidence["whole_reference_verified"] is False
        assert client.post("/api/sanger-verifications", json={"sequencing_read_id": uploaded.json()["id"]}).status_code == 409
    with TestClient(app) as client:
        saved = client.get(f"/api/sequence-revisions/{revision}/sanger-reads").json()[0]
        assert saved["evidence"] == evidence
        assert "storage_path" not in saved
        assert client.get("/api/sequence-revisions/missing/sanger-reads").status_code == 404
        read_id = uploaded.json()["id"]
        trace = client.get(f"/api/sanger-reads/{read_id}/trace?start=2&count=8")
        assert trace.status_code == 200
        assert trace.json()["available"] is True
        assert len(trace.json()["bases"]) == 8
        assert trace.json()["bases"][0]["base"] == reference[-28]
        assert client.get(f"/api/sanger-reads/{read_id}/trace?start=-1").status_code == 422
        assert client.get("/api/sanger-reads/missing/trace").status_code == 404
        with connect() as db:
            job = db.execute("SELECT * FROM analysis_jobs WHERE job_kind = 'sanger-verification'").fetchone()
            assert json.loads(job["parameters_json"])["evidence_version"] == 4
            assert json.loads(job["input_manifest_json"])["read_sha256"] == uploaded.json()["file_sha256"]
            assert len(json.loads(db.execute("SELECT qualities_json FROM sequencing_reads").fetchone()[0])) == 30
            stored = db.execute("SELECT storage_path FROM sequencing_reads").fetchone()[0]
        from pathlib import Path
        Path(stored).write_bytes(b"changed fixture")
        assert client.get(f"/api/sanger-reads/{read_id}/trace").status_code == 409
        Path(stored).unlink()
        assert client.get(f"/api/sanger-reads/{read_id}/trace").status_code == 409


def test_missing_revision_rejected_before_ab1_parse(tmp_path, monkeypatch):
    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    def unexpected_parse(_):
        raise AssertionError("Should validate revision before parsing or storing read")
    monkeypatch.setattr("localmolbio.api.parse_ab1", unexpected_parse)
    with TestClient(app) as client:
        response = client.post("/api/sanger-reads", data={"sequence_revision_id": "missing"}, files={"file": ("trace.ab1", b"fixture")})
        assert response.status_code == 404
        assert not (tmp_path / "runtime" / "reads").exists()


def test_repaired_origin_annotations_keep_raw_text_both_strands_and_warning_details(tmp_path, monkeypatch):
    import json
    from io import StringIO
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    reference = SeqRecord(Seq('ACGT' * 3), id='wrapped', annotations={'molecule_type':'DNA','topology':'circular'})
    reference.features = [
        SeqFeature(CompoundLocation([FeatureLocation(9,12,strand=1),FeatureLocation(0,3,strand=1)]),type='misc_feature'),
        SeqFeature(CompoundLocation([FeatureLocation(0,2,strand=-1),FeatureLocation(8,12,strand=-1)]),type='misc_feature')]
    stream = StringIO();SeqIO.write(reference, stream, 'genbank')
    record = stream.getvalue().replace('join(10..12,1..3)', '10..3').replace('complement(join(9..12,1..2))','complement(9..2)')
    archive = BytesIO()
    with ZipFile(archive, 'w') as z:
        z.writestr('synthetic.gb', record)
    monkeypatch.setenv('MOLBIO_DATA_DIR', str(tmp_path / 'runtime'))
    with TestClient(app) as client:
        imported = client.post('/api/imports/benchling', files={'file': ('synthetic.zip', archive.getvalue(), 'application/zip')})
        assert imported.status_code == 201
        assert imported.json()['parse_warning_count'] == 2
        sequence = client.get('/api/sequences').json()[0]['id']
        expected = [[{'start':9,'end':12,'strand':1},{'start':0,'end':3,'strand':1}],
                    [{'start':0,'end':2,'strand':-1},{'start':8,'end':12,'strand':-1}]]
        for legacy in [False, True]:
            if legacy:
                with connect() as db: db.execute('UPDATE sequences SET features_json = ?', ('[]',))
            result = client.get(f'/api/sequences/{sequence}/map').json()
            assert result['requires_annotation_review']
            assert len(json.loads(result['parse_warnings_json'])) == 2
            assert [feature['segments'] for feature in result['features']] == expected
            with connect() as db:
                assert db.execute('SELECT raw_genbank FROM sequences WHERE id=?', (sequence,)).fetchone()[0] == record
