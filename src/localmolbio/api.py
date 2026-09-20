from __future__ import annotations

from contextlib import asynccontextmanager
import csv
from datetime import datetime, timezone
from io import StringIO
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal, Optional
from uuid import uuid4

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel, Field

from .config import database_path
from .database import connect, initialise
from .importer import (
    ImportErrorDetail,
    _parse_genbank,
    ensure_initial_revisions,
    import_benchling_archive,
)
from .primer_design import PrimerDesignError, PrimerDesignSettings, design_pcr_primers


class PrimerDesignInput(BaseModel):
    sequence_revision_id: str
    purpose: Literal["pcr"] = "pcr"
    name_prefix: Optional[str] = Field(default=None, max_length=100)
    product_size_min: int = Field(default=100, ge=50, le=10_000)
    product_size_max: int = Field(default=800, ge=51, le=10_000)
    num_return: int = Field(default=5, ge=1, le=25)
    min_tm: float = Field(default=57.0, ge=40.0, le=80.0)
    opt_tm: float = Field(default=60.0, ge=40.0, le=80.0)
    max_tm: float = Field(default=63.0, ge=40.0, le=80.0)


class PrimerSelectionInput(BaseModel):
    selection_state: Literal["selected", "archived"]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@asynccontextmanager
async def lifespan(_: FastAPI):
    initialise()
    ensure_initial_revisions()
    yield


app = FastAPI(
    title="Local Molecular Biology Workbench", version="0.1.0", lifespan=lifespan
)
STATIC_DIR = Path(__file__).parent / "static"


@app.get("/", response_class=FileResponse)
def index() -> Path:
    return STATIC_DIR / "index.html"


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/api/imports")
def list_imports() -> list[dict[str, object]]:
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT id, original_filename, archive_sha256, imported_at, total_members,
                   imported_records, skipped_members, parse_warning_records,
                   parse_warning_count
            FROM imports ORDER BY imported_at DESC
            """
        ).fetchall()
    return [dict(row) for row in rows]


@app.get("/api/dashboard")
def dashboard() -> dict[str, int]:
    with connect() as connection:
        row = connection.execute(
            """
            SELECT
                COUNT(*) AS sequence_count,
                COALESCE(SUM(topology = 'circular'), 0) AS circular_count,
                COALESCE(SUM(parse_warning_count > 0), 0) AS review_count
            FROM sequences
            """
        ).fetchone()
        selected_primers = connection.execute(
            "SELECT COUNT(*) AS count FROM primers WHERE selection_state = 'selected'"
        ).fetchone()["count"]
    return {**dict(row), "selected_primer_count": selected_primers}


@app.get("/api/sequences")
def list_sequences(query: str = "", limit: int = 200) -> list[dict[str, object]]:
    safe_limit = min(max(limit, 1), 500)
    like = f"%{query.strip()}%"
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT id, display_name, archive_member_name, archive_member_occurrence,
                   length_bp, topology, molecule_type, feature_count, sequence_sha256
                   , parse_warning_count
            FROM sequences
            WHERE display_name LIKE ? OR archive_member_name LIKE ?
            ORDER BY display_name COLLATE NOCASE, archive_member_index
            LIMIT ?
            """,
            (like, like, safe_limit),
        ).fetchall()
    return [dict(row) for row in rows]


@app.get("/api/sequences/{sequence_id}")
def sequence_detail(sequence_id: str) -> dict[str, object]:
    with connect() as connection:
        row = connection.execute(
            """
            SELECT id, import_id, archive_member_index, archive_member_name,
                   archive_member_occurrence, member_sha256, sequence_sha256,
                   display_name, length_bp, topology, molecule_type, feature_count,
                   features_json, parse_warning_count, parse_warnings_json, created_at
            FROM sequences WHERE id = ?
            """,
            (sequence_id,),
        ).fetchone()
    if not row:
        raise HTTPException(status_code=404, detail="Sequence not found")
    return dict(row)


@app.get("/api/sequences/{sequence_id}/revisions")
def list_sequence_revisions(sequence_id: str) -> list[dict[str, object]]:
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT id, parent_revision_id, revision_number, label, sequence_sha256,
                   topology, source_kind, created_at
            FROM sequence_revisions
            WHERE sequence_id = ?
            ORDER BY revision_number DESC
            """,
            (sequence_id,),
        ).fetchall()
    if not rows:
        raise HTTPException(status_code=404, detail="Sequence not found")
    return [dict(row) for row in rows]


@app.get("/api/jobs")
def list_analysis_jobs(limit: int = 100) -> list[dict[str, object]]:
    safe_limit = min(max(limit, 1), 500)
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT id, sequence_revision_id, job_kind, status, created_at, started_at,
                   completed_at, error_detail
            FROM analysis_jobs
            ORDER BY created_at DESC LIMIT ?
            """,
            (safe_limit,),
        ).fetchall()
    return [dict(row) for row in rows]


@app.get("/api/revisions/{revision_id}/primers")
def list_primers(revision_id: str) -> list[dict[str, object]]:
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT id, name, sequence_text, direction, purpose, binding_start,
                   binding_end, metrics_json, design_parameters_json,
                   selection_state, created_at
            FROM primers
            WHERE sequence_revision_id = ?
            ORDER BY created_at DESC, name COLLATE NOCASE
            """,
            (revision_id,),
        ).fetchall()
        revision_exists = connection.execute(
            "SELECT 1 FROM sequence_revisions WHERE id = ?", (revision_id,)
        ).fetchone()
    if not revision_exists:
        raise HTTPException(status_code=404, detail="Sequence revision not found")
    primers = []
    for row in rows:
        item = dict(row)
        item["metrics"] = json.loads(item.pop("metrics_json"))
        item["design_parameters"] = json.loads(item.pop("design_parameters_json"))
        primers.append(item)
    return primers


@app.patch("/api/primers/{primer_id}")
def update_primer_selection(
    primer_id: str, update: PrimerSelectionInput
) -> dict[str, object]:
    now = _utc_now()
    with connect() as connection:
        primer = connection.execute(
            """
            SELECT id, sequence_revision_id, name, selection_state
            FROM primers WHERE id = ?
            """,
            (primer_id,),
        ).fetchone()
        if not primer:
            raise HTTPException(status_code=404, detail="Primer not found")
        connection.execute(
            "UPDATE primers SET selection_state = ? WHERE id = ?",
            (update.selection_state, primer_id),
        )
        connection.execute(
            """
            INSERT INTO audit_events (id, object_type, object_id, action, payload_json, created_at)
            VALUES (?, 'primer', ?, 'primer.selection-updated', ?, ?)
            """,
            (
                str(uuid4()),
                primer_id,
                json.dumps(
                    {
                        "from": primer["selection_state"],
                        "to": update.selection_state,
                        "sequence_revision_id": primer["sequence_revision_id"],
                    }
                ),
                now,
            ),
        )
    return {"id": primer_id, "selection_state": update.selection_state}


@app.get("/api/revisions/{revision_id}/primers.csv")
def export_selected_primers(revision_id: str) -> Response:
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT name, sequence_text, direction, purpose, binding_start, binding_end,
                   metrics_json, created_at
            FROM primers
            WHERE sequence_revision_id = ? AND selection_state = 'selected'
            ORDER BY name COLLATE NOCASE
            """,
            (revision_id,),
        ).fetchall()
        revision_exists = connection.execute(
            "SELECT 1 FROM sequence_revisions WHERE id = ?", (revision_id,)
        ).fetchone()
    if not revision_exists:
        raise HTTPException(status_code=404, detail="Sequence revision not found")

    output = StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=[
            "name",
            "sequence",
            "direction",
            "purpose",
            "binding_start_1_based",
            "binding_end_1_based_inclusive",
            "tm_celsius",
            "gc_percent",
            "source_revision_id",
            "created_at",
        ],
    )
    writer.writeheader()
    for row in rows:
        metrics = json.loads(row["metrics_json"])
        writer.writerow(
            {
                "name": row["name"] or "",
                "sequence": row["sequence_text"],
                "direction": row["direction"],
                "purpose": row["purpose"],
                "binding_start_1_based": row["binding_start"] + 1,
                "binding_end_1_based_inclusive": row["binding_end"],
                "tm_celsius": metrics.get("tm", ""),
                "gc_percent": metrics.get("gc_percent", ""),
                "source_revision_id": revision_id,
                "created_at": row["created_at"],
            }
        )
    return Response(
        content=output.getvalue(),
        media_type="text/csv; charset=utf-8",
        headers={
            "Content-Disposition": f'attachment; filename="selected-primers-{revision_id}.csv"'
        },
    )


@app.post("/api/primer-designs", status_code=201)
def create_primer_design(request: PrimerDesignInput) -> dict[str, object]:
    settings = PrimerDesignSettings(
        product_size_min=request.product_size_min,
        product_size_max=request.product_size_max,
        num_return=request.num_return,
        min_tm=request.min_tm,
        opt_tm=request.opt_tm,
        max_tm=request.max_tm,
    )
    parameters = request.model_dump()
    job_id = str(uuid4())
    now = _utc_now()
    with connect() as connection:
        revision = connection.execute(
            """
            SELECT id, label, sequence_text, sequence_sha256
            FROM sequence_revisions WHERE id = ?
            """,
            (request.sequence_revision_id,),
        ).fetchone()
        if not revision:
            raise HTTPException(status_code=404, detail="Sequence revision not found")
        connection.execute(
            """
            INSERT INTO analysis_jobs (
                id, sequence_revision_id, job_kind, status, parameters_json,
                input_manifest_json, result_summary_json, created_at, started_at
            ) VALUES (?, ?, 'primer-design', 'running', ?, ?, '{}', ?, ?)
            """,
            (
                job_id,
                revision["id"],
                json.dumps(parameters),
                json.dumps({"sequence_sha256": revision["sequence_sha256"]}),
                now,
                now,
            ),
        )
    try:
        pairs = design_pcr_primers(revision["sequence_text"], settings)
    except PrimerDesignError as exc:
        with connect() as connection:
            connection.execute(
                """
                UPDATE analysis_jobs
                SET status = 'failed', error_detail = ?, completed_at = ?
                WHERE id = ?
                """,
                (str(exc), _utc_now(), job_id),
            )
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    created_primers: list[dict[str, object]] = []
    prefix = request.name_prefix or revision["label"]
    with connect() as connection:
        for pair in pairs:
            for direction, item, suffix in (
                ("forward", pair["left"], "F"),
                ("reverse", pair["right"], "R"),
            ):
                primer_id = str(uuid4())
                metrics = {
                    "tm": item["tm"],
                    "gc_percent": item["gc_percent"],
                    "self_any_th": item["self_any_th"],
                    "self_end_th": item["self_end_th"],
                    "product_size": pair["product_size"],
                    "pair_index": pair["pair_index"],
                }
                name = f"{prefix}-{suffix}{pair['pair_index']}"
                connection.execute(
                    """
                    INSERT INTO primers (
                        id, sequence_revision_id, name, sequence_text, direction,
                        purpose, binding_start, binding_end, metrics_json,
                        design_parameters_json, selection_state, created_at
                    ) VALUES (?, ?, ?, ?, ?, 'pcr', ?, ?, ?, ?, 'candidate', ?)
                    """,
                    (
                        primer_id,
                        revision["id"],
                        name,
                        item["sequence"],
                        direction,
                        item["binding_start"],
                        item["binding_end"],
                        json.dumps(metrics),
                        json.dumps(parameters),
                        now,
                    ),
                )
                created_primers.append(
                    {
                        "id": primer_id,
                        "name": name,
                        "direction": direction,
                        "sequence": item["sequence"],
                        "binding_start": item["binding_start"],
                        "binding_end": item["binding_end"],
                        "metrics": metrics,
                    }
                )
        summary = {"pair_count": len(pairs), "primer_count": len(created_primers)}
        connection.execute(
            """
            UPDATE analysis_jobs
            SET status = 'succeeded', result_summary_json = ?, completed_at = ?
            WHERE id = ?
            """,
            (json.dumps(summary), _utc_now(), job_id),
        )
        connection.execute(
            """
            INSERT INTO audit_events (id, object_type, object_id, action, payload_json, created_at)
            VALUES (?, 'analysis_job', ?, 'primer-design.completed', ?, ?)
            """,
            (str(uuid4()), job_id, json.dumps(summary), _utc_now()),
        )
    return {"job_id": job_id, "pairs": pairs, "primers": created_primers}


@app.get("/api/sequences/{sequence_id}/map")
def sequence_map(sequence_id: str) -> dict[str, object]:
    with connect() as connection:
        row = connection.execute(
            """
            SELECT id, display_name, length_bp, topology, feature_count,
                   features_json, parse_warning_count, parse_warnings_json,
                   archive_member_name, raw_genbank,
                   (
                       SELECT id FROM sequence_revisions
                       WHERE sequence_id = sequences.id
                       ORDER BY revision_number DESC LIMIT 1
                   ) AS current_revision_id
            FROM sequences WHERE id = ?
            """,
            (sequence_id,),
        ).fetchone()
    if not row:
        raise HTTPException(status_code=404, detail="Sequence not found")

    result = dict(row)
    features = json.loads(result.pop("features_json"))
    if not features and result["feature_count"]:
        features = _parse_genbank(
            result.pop("raw_genbank"), result["archive_member_name"]
        )["features"]
    else:
        result.pop("raw_genbank")
    result["features"] = features
    result["requires_annotation_review"] = bool(result["parse_warning_count"])
    return result


@app.post("/api/imports/benchling", status_code=201)
async def import_benchling_export(file: UploadFile = File(...)) -> dict[str, object]:
    if not file.filename or not file.filename.lower().endswith(".zip"):
        raise HTTPException(status_code=400, detail="Upload a Benchling ZIP export")

    with NamedTemporaryFile(suffix=".zip", delete=False) as temporary:
        temporary_path = Path(temporary.name)
        while block := await file.read(1024 * 1024):
            temporary.write(block)

    try:
        result = import_benchling_archive(temporary_path, database_path())
    except ImportErrorDetail as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    finally:
        temporary_path.unlink(missing_ok=True)

    return {
        "import_id": result.import_id,
        "total_members": result.total_members,
        "imported_records": result.imported_records,
        "skipped_members": result.skipped_members,
        "parse_warning_records": result.parse_warning_records,
        "parse_warning_count": result.parse_warning_count,
        "archive_sha256": result.archive_sha256,
    }
