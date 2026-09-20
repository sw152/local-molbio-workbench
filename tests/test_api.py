from __future__ import annotations

from io import BytesIO
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi.testclient import TestClient

from localmolbio.api import app


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

        detail = client.get(f"/api/sequences/{sequences.json()[0]['id']}")
        assert detail.status_code == 200
        assert "spans the origin" in detail.json()["parse_warnings_json"]
