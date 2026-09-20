from __future__ import annotations

from io import BytesIO
import random
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi.testclient import TestClient

from localmolbio.api import app


def test_pcr_primer_design_persists_candidates_and_job(tmp_path, monkeypatch) -> None:
    random_source = random.Random(42)
    bases = "".join(random_source.choice("ACGT") for _ in range(1800))
    record = f"""LOCUS       primer_test            1800 bp    ds-DNA     linear       01-JAN-2026
DEFINITION  synthetic primer-design test record.
FEATURES             Location/Qualifiers
     source          1..1800
                     /organism=\"synthetic construct\"
ORIGIN
        1 {bases.lower()}
//
"""
    archive = BytesIO()
    with ZipFile(archive, "w", ZIP_DEFLATED) as zipper:
        zipper.writestr("primer-test.gb", record)

    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    with TestClient(app) as client:
        imported = client.post(
            "/api/imports/benchling",
            files={"file": ("benchling.zip", archive.getvalue(), "application/zip")},
        )
        assert imported.status_code == 201
        sequence_id = client.get("/api/sequences").json()[0]["id"]
        revision_id = client.get(f"/api/sequences/{sequence_id}/revisions").json()[0]["id"]

        designed = client.post(
            "/api/primer-designs",
            json={
                "sequence_revision_id": revision_id,
                "name_prefix": "demo",
                "product_size_min": 150,
                "product_size_max": 300,
                "num_return": 3,
            },
        )
        assert designed.status_code == 201
        body = designed.json()
        assert len(body["pairs"]) >= 1
        assert len(body["primers"]) == len(body["pairs"]) * 2
        assert body["primers"][0]["name"] == "demo-F1"

        saved_primers = client.get(f"/api/revisions/{revision_id}/primers")
        assert saved_primers.status_code == 200
        assert len(saved_primers.json()) == len(body["primers"])
        assert saved_primers.json()[0]["metrics"]["tm"] > 40

        selected = client.patch(
            f"/api/primers/{body['primers'][0]['id']}",
            json={"selection_state": "selected"},
        )
        assert selected.status_code == 200
        assert selected.json()["selection_state"] == "selected"
        exported = client.get(f"/api/revisions/{revision_id}/primers.csv")
        assert exported.status_code == 200
        assert "demo-F1" in exported.text
        assert "source_revision_id" in exported.text

        jobs = client.get("/api/jobs").json()
        assert jobs[0]["id"] == body["job_id"]
        assert jobs[0]["status"] == "succeeded"
