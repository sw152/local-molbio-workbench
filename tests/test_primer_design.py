from __future__ import annotations

from io import BytesIO, StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi.testclient import TestClient

from localmolbio.api import app


def test_pcr_primer_design_persists_candidates_and_job(tmp_path, monkeypatch) -> None:
    random_source = random.Random(42)
    bases = "".join(random_source.choice("ACGT") for _ in range(1800))
    record_object = SeqRecord(Seq(bases), id="primer_test", name="primer_test",
                              annotations={"molecule_type":"DNA", "topology":"linear"})
    stream = StringIO()
    SeqIO.write(record_object, stream, "genbank")
    record = stream.getvalue()
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

        targeted = client.post("/api/primer-designs", json={
            "sequence_revision_id": revision_id, "product_size_min": 300, "product_size_max": 500,
            "num_return": 2, "target_start": 650, "target_end": 850,
        })
        assert targeted.status_code == 201, targeted.text
        evidence = targeted.json()
        assert evidence["pairs"]
        assert evidence["parameters"]["reference_topology"] == "linear"
        assert evidence["parameters"]["target_start"] == 650
        assert evidence["parameters"]["specificity_status"] == "not_evaluated"
        assert evidence["parameters"]["engine_version"]
        saved = client.get(f"/api/revisions/{revision_id}/primers").json()
        matching = [p for p in saved if p["metrics"].get("analysis_job_id") == evidence["job_id"]]
        assert len(matching) == len(evidence["primers"])
        assert matching[0]["metrics"]["target"] == {"start":650,"end":850}
        assert matching[0]["design_parameters"]["target_end"] == 850
        import csv
        client.patch(f"/api/primers/{matching[0]['id']}",json={"selection_state":"selected"})
        rows = list(csv.DictReader(StringIO(client.get(f"/api/revisions/{revision_id}/primers.csv").text)))
        target_row = next(r for r in rows if r["analysis_job_id"] == evidence["job_id"])
        assert target_row["target_start_1_based"] == "651"
        assert target_row["target_end_1_based_inclusive"] == "850"
        assert target_row["specificity_status"] == "not_evaluated"
        invalid = client.post("/api/primer-designs",json={"sequence_revision_id":revision_id,"target_start":650})
        assert invalid.status_code == 422
        assert client.get("/api/jobs").json()[0]["status"] == "failed"
        def broken_engine(*args):
            raise RuntimeError("simulated engine failure")
        monkeypatch.setattr("localmolbio.api.design_pcr_primers",broken_engine)
        failed = client.post("/api/primer-designs",json={"sequence_revision_id":revision_id})
        assert failed.status_code == 500
        assert all(job["status"] != "running" for job in client.get("/api/jobs").json())
        assert len(client.get(f"/api/revisions/{revision_id}/primers").json()) == len(saved)



import pytest
from Bio.Seq import Seq
from localmolbio.primer_design import PrimerDesignError, PrimerDesignSettings, design_pcr_primers


def synthetic_template():
    rng=random.Random(42)
    return ''.join(rng.choice('ACGT') for _ in range(1800))


@pytest.mark.parametrize('target',[(650,850),(700,701)])
def test_target_pairs_flank_interval_and_oligos_match_reference(target):
    template=synthetic_template()
    pairs=design_pcr_primers(template,PrimerDesignSettings(product_size_min=300,product_size_max=500,target_start=target[0],target_end=target[1]))
    assert pairs
    for pair in pairs:
        left,right=pair['left'],pair['right']
        assert left['binding_end']<=target[0]<target[1]<=right['binding_start']
        assert template[left['binding_start']:left['binding_end']]==left['sequence']
        assert str(Seq(template[right['binding_start']:right['binding_end']]).reverse_complement())==right['sequence']
        assert right['binding_end']-left['binding_start']==pair['product_size']
        assert pair['specificity_status']=='not_evaluated'


@pytest.mark.parametrize('options',[
    {'target_start':100}, {'target_end':200}, {'target_start':200,'target_end':100},
    {'target_start':0,'target_end':50}, {'target_start':1750,'target_end':1800},
    {'target_start':200,'target_end':200}, {'target_start':100,'target_end':1801},
    {'target_start':100,'target_end':900,'product_size_max':800},
    {'min_tm':64,'opt_tm':60}, {'opt_tm':70,'max_tm':63},
    {'min_tm':float('nan')}, {'num_return':0}, {'product_size_min':49},
])
def test_invalid_target_and_parameter_constraints(options):
    with pytest.raises(PrimerDesignError):
        design_pcr_primers(synthetic_template(),PrimerDesignSettings(**options))


def test_no_candidates_is_distinct_from_invalid_target():
    assert design_pcr_primers('A'*1000,PrimerDesignSettings(target_start=300,target_end=400))==[]
