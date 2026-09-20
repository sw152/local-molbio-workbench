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
