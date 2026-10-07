from __future__ import annotations

from contextlib import asynccontextmanager
import csv
from datetime import datetime, timezone
from io import StringIO
import json
from pathlib import Path
from shutil import copyfile
from tempfile import NamedTemporaryFile
from typing import Literal, Optional
from uuid import uuid4

from fastapi import FastAPI, File, Form, HTTPException, UploadFile, Query
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .config import database_path, read_dir
from .chromatogram import trace_window
from .database import connect, initialise
from .importer import (
    ImportErrorDetail,
    _parse_genbank,
    ensure_initial_revisions,
    import_benchling_archive,
)
from .primer_sites import exact_reference_sites
from .primer_design import PrimerDesignError, PrimerDesignSettings, design_pcr_primers
import primer3
from .sanger import SangerReadError, file_sha256, parse_ab1
from .sanger_verification import ALIGNMENT_PARAMETERS, SangerVerificationError, align_sanger_read


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
    target_start: Optional[int] = Field(default=None, ge=0)
    target_end: Optional[int] = Field(default=None, ge=1)


class PrimerSelectionInput(BaseModel):
    selection_state: Literal["selected", "archived"]


class SangerVerificationInput(BaseModel):
    sequencing_read_id: str
    previous_alignment_id: Optional[str] = None
    direction: Optional[Literal["forward", "reverse", "unknown"]] = None


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
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


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


@app.get("/api/sequence-library")
def paginated_sequences(
    query: str = Query(default="", max_length=500),
    limit: int = Query(default=24, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
) -> dict[str, object]:
    # Escape LIKE metacharacters: a name containing % or _ is a literal search.
    escaped = query.strip().replace("!", "!!").replace("%", "!%").replace("_", "!_")
    like = f"%{escaped}%"
    where = "display_name LIKE ? ESCAPE '!' OR archive_member_name LIKE ? ESCAPE '!'"
    with connect() as connection:
        # The count and page share a read snapshot even if another client imports.
        connection.execute("BEGIN")
        total = connection.execute(f"SELECT COUNT(*) FROM sequences WHERE {where}", (like, like)).fetchone()[0]
        rows = connection.execute(
            f"""SELECT id, display_name, archive_member_name, archive_member_occurrence,
                       length_bp, topology, molecule_type, feature_count, sequence_sha256, parse_warning_count
                FROM sequences WHERE {where}
                ORDER BY display_name COLLATE NOCASE, archive_member_index, id
                LIMIT ? OFFSET ?""", (like, like, limit, offset),
        ).fetchall()
    return {"items": [dict(row) for row in rows], "total": total, "offset": offset,
            "limit": limit, "has_more": offset + len(rows) < total}


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
            "analysis_job_id",
            "pair_index",
            "product_size_bp",
            "target_start_1_based",
            "target_end_1_based_inclusive",
            "specificity_status",
            "reference_exact_directional_matches",
            "reference_site_review_method",
            "created_at",
        ],
    )
    writer.writeheader()
    for row in rows:
        metrics = json.loads(row["metrics_json"])
        target = metrics.get("target") or {}
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
                "analysis_job_id": metrics.get("analysis_job_id", ""),
                "pair_index": metrics.get("pair_index", ""),
                "product_size_bp": metrics.get("product_size", ""),
                "target_start_1_based": target["start"] + 1 if "start" in target else "",
                "target_end_1_based_inclusive": target.get("end", ""),
                "specificity_status": metrics.get("specificity_status", "not_recorded"),
                "reference_exact_directional_matches": metrics.get("reference_sites", {}).get("total_directional_matches", ""),
                "reference_site_review_method": metrics.get("reference_sites", {}).get("method", "not_recorded"),
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
        target_start=request.target_start,
        target_end=request.target_end,
    )
    parameters = {**request.model_dump(), "engine": "primer3-py", "engine_version": primer3.__version__,
                  "coordinate_system": "zero-based-half-open", "design_scope": "linear_template_coordinates",
                  "specificity_status": "not_evaluated", "origin_spanning_supported": False}
    job_id = str(uuid4())
    now = _utc_now()
    with connect() as connection:
        revision = connection.execute(
            """
            SELECT id, label, sequence_text, sequence_sha256, topology
            FROM sequence_revisions WHERE id = ?
            """,
            (request.sequence_revision_id,),
        ).fetchone()
        if not revision:
            raise HTTPException(status_code=404, detail="Sequence revision not found")
        parameters["reference_topology"] = revision["topology"]
        parameters["binding_site_review"] = "full_length_exact_match_v1"
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
        for pair in pairs:
            for side in ("left", "right"):
                pair[side]["reference_sites"] = exact_reference_sites(
                    revision["sequence_text"], pair[side]["sequence"], revision["topology"])
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

    except Exception as exc:
        with connect() as connection:
            connection.execute(
                "UPDATE analysis_jobs SET status = 'failed', error_detail = ?, completed_at = ? WHERE id = ?",
                (f"{type(exc).__name__}: {exc}", _utc_now(), job_id),
            )
        raise HTTPException(status_code=500, detail="Primer design engine failed; no candidates were saved") from exc

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
                    "reference_sites": item["reference_sites"],
                    "product_size": pair["product_size"],
                    "pair_index": pair["pair_index"],
                    "analysis_job_id": job_id,
                    "target": pair["target"],
                    "specificity_status": "not_evaluated",
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
        summary = {"pair_count": len(pairs), "primer_count": len(created_primers), "parameters": parameters}
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
    return {"job_id": job_id, "pairs": pairs, "primers": created_primers, "parameters": parameters}


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


@app.post("/api/sanger-reads", status_code=201)
async def upload_sanger_read(
    sequence_revision_id: str = Form(...),
    direction: Literal["forward", "reverse", "unknown"] = Form("unknown"),
    file: UploadFile = File(...),
) -> dict[str, object]:
    if not file.filename or not file.filename.lower().endswith(".ab1"):
        raise HTTPException(status_code=400, detail="Upload an AB1 chromatogram file")
    with connect() as connection:
        if not connection.execute("SELECT id FROM sequence_revisions WHERE id = ?", (sequence_revision_id,)).fetchone():
            raise HTTPException(status_code=404, detail="Sequence revision not found")
    with NamedTemporaryFile(suffix=".ab1", delete=False) as temporary:
        temporary_path = Path(temporary.name)
        while block := await file.read(1024 * 1024):
            temporary.write(block)
    try:
        read_hash = file_sha256(temporary_path)
        parsed = parse_ab1(temporary_path)
        destination_directory = read_dir()
        destination_directory.mkdir(parents=True, exist_ok=True)
        destination = destination_directory / f"{read_hash}.ab1"
        if not destination.exists():
            copyfile(temporary_path, destination)
        read_id = str(uuid4())
        now = _utc_now()
        with connect() as connection:
            revision = connection.execute(
                "SELECT id FROM sequence_revisions WHERE id = ?", (sequence_revision_id,)
            ).fetchone()
            if not revision:
                raise HTTPException(status_code=404, detail="Sequence revision not found")
            duplicate = connection.execute(
                "SELECT id FROM sequencing_reads WHERE sequence_revision_id = ? AND file_sha256 = ?",
                (sequence_revision_id, read_hash),
            ).fetchone()
            if duplicate:
                raise HTTPException(status_code=409, detail="This AB1 read is already attached to the revision")
            connection.execute(
                """
                INSERT INTO sequencing_reads (
                    id, sequence_revision_id, original_filename, file_sha256, storage_path,
                    file_format, direction, base_sequence, length_bp, quality_summary_json,
                    parser_metadata_json, created_at, qualities_json
                ) VALUES (?, ?, ?, ?, ?, 'ab1', ?, ?, ?, ?, ?, ?, ?)
                """,
                (read_id, sequence_revision_id, file.filename, read_hash, str(destination), direction,
                 parsed.sequence, len(parsed.sequence), json.dumps(parsed.quality_summary),
                 json.dumps(parsed.parser_metadata), now, json.dumps(parsed.qualities)),
            )
        return {"id": read_id, "length_bp": len(parsed.sequence), "direction": direction,
                "quality_summary": parsed.quality_summary, "file_sha256": read_hash}
    except SangerReadError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        temporary_path.unlink(missing_ok=True)


def _analysis_report(row) -> dict:
    report = json.loads(row["report_json"])
    if not report.get("evidence"):
        report["evidence"] = {"scope": "legacy_local_alignment", "review_flags": ["legacy_evidence_unavailable"], "whole_reference_verified": False}
    return {**report, "alignment_id": row["id"], "job_id": row["analysis_job_id"],
            "run_number": row["run_number"], "previous_alignment_id": row["previous_alignment_id"],
            "analyzed_at": row["created_at"]}


def _check_previous(connection, read_id, expected):
    latest = connection.execute("SELECT id, run_number FROM sanger_analysis_runs WHERE sequencing_read_id = ? ORDER BY run_number DESC LIMIT 1", (read_id,)).fetchone()
    if (latest["id"] if latest else None) != expected:
        raise HTTPException(status_code=409, detail="Analysis history changed or a report already exists; reload before starting a new analysis")
    return latest["run_number"] + 1 if latest else 1


@app.post("/api/sanger-verifications", status_code=201)
def create_sanger_verification(request: SangerVerificationInput) -> dict[str, object]:
    with connect() as connection:
        row = connection.execute(
            """SELECT reads.*, revisions.sequence_text, revisions.sequence_sha256, revisions.topology
               FROM sequencing_reads AS reads
               JOIN sequence_revisions AS revisions ON revisions.id = reads.sequence_revision_id
               WHERE reads.id = ?""", (request.sequencing_read_id,),
        ).fetchone()
        if not row:
            raise HTTPException(status_code=404, detail="Sanger read not found")
        _check_previous(connection, row["id"], request.previous_alignment_id)
    # Reparse hash-checked original input; historical reads may lack saved per-base qualities.
    try:
        path = Path(row["storage_path"])
        if file_sha256(path) != row["file_sha256"]:
            raise HTTPException(status_code=409, detail="Original AB1 hash changed; analysis was not saved")
        parsed = parse_ab1(path)
        if file_sha256(path) != row["file_sha256"]:
            raise HTTPException(status_code=409, detail="Original AB1 changed during parsing; analysis was not saved")
    except OSError as exc:
        raise HTTPException(status_code=409, detail="Original AB1 file is unavailable; analysis was not saved") from exc
    except SangerReadError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    saved_qualities = json.loads(row["qualities_json"])
    if parsed.sequence != row["base_sequence"] or (saved_qualities and saved_qualities != parsed.qualities):
        raise HTTPException(status_code=409, detail="Original AB1 calls or qualities disagree with the saved read; analysis was not saved")
    direction = request.direction or row["direction"]
    started = _utc_now()
    try:
        result = align_sanger_read(row["sequence_text"], parsed.sequence, direction,
                                   row["topology"] == "circular", qualities=parsed.qualities or None)
    except SangerVerificationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail="Alignment engine failed; existing reports were preserved") from exc
    job_id, alignment_id, now = str(uuid4()), str(uuid4()), _utc_now()
    report = {key: getattr(result, key) for key in ("reference_start", "reference_end", "wraps_origin", "aligned_bases", "matched_bases", "mismatched_bases", "inserted_bases", "deleted_bases", "identity_fraction", "variants", "evidence")}
    parameters = {**ALIGNMENT_PARAMETERS, "requested_direction": direction,
                  "original_read_direction": row["direction"], "quality_source": "hash_checked_original_ab1",
                  "reference_topology": row["topology"], "trimming": "none"}
    result.evidence["parameters"] = parameters
    summary = {"identity_fraction": result.identity_fraction, "variant_count": len(result.variants), "evidence": result.evidence}
    with connect() as connection:
        # Check again under the write lock: simultaneous requests cannot create two runs.
        connection.execute("BEGIN IMMEDIATE")
        run_number = _check_previous(connection, row["id"], request.previous_alignment_id)
        connection.execute(
            """INSERT INTO analysis_jobs (id, sequence_revision_id, job_kind, status, parameters_json,
                input_manifest_json, result_summary_json, created_at, started_at, completed_at)
               VALUES (?, ?, 'sanger-verification', 'succeeded', ?, ?, ?, ?, ?, ?)""",
            (job_id, row["sequence_revision_id"], json.dumps(parameters),
             json.dumps({"read_id": row["id"], "read_sha256": row["file_sha256"], "reference_sha256": row["sequence_sha256"], "previous_alignment_id": request.previous_alignment_id}),
             json.dumps(summary), started, started, now),
        )
        connection.execute(
            """INSERT INTO sanger_analysis_runs (id, sequencing_read_id, sequence_revision_id, analysis_job_id,
               run_number, previous_alignment_id, report_json, created_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (alignment_id, row["id"], row["sequence_revision_id"], job_id, run_number, request.previous_alignment_id, json.dumps(report), now),
        )
        connection.execute("INSERT INTO audit_events (id, object_type, object_id, action, payload_json, created_at) VALUES (?, 'sequencing-read', ?, 'analysis-created', ?, ?)",
                           (str(uuid4()), row["id"], json.dumps({"alignment_id":alignment_id,"previous_alignment_id":request.previous_alignment_id,"run_number":run_number}), now))
    return {"id": alignment_id, "job_id": job_id, "run_number": run_number,
            "previous_alignment_id": request.previous_alignment_id, **summary, **report}


@app.get("/api/sequence-revisions/{revision_id}/sanger-reads")
def list_sanger_reads(revision_id: str) -> list[dict[str, object]]:
    with connect() as connection:
        if not connection.execute("SELECT id FROM sequence_revisions WHERE id = ?", (revision_id,)).fetchone():
            raise HTTPException(status_code=404, detail="Sequence revision not found")
        rows = connection.execute("SELECT id, original_filename, length_bp, direction, quality_summary_json, created_at FROM sequencing_reads WHERE sequence_revision_id = ? ORDER BY created_at DESC", (revision_id,)).fetchall()
        output = []
        for row in rows:
            item = dict(row)
            item["quality_summary"] = json.loads(item.pop("quality_summary_json"))
            latest = connection.execute("SELECT * FROM sanger_analysis_runs WHERE sequencing_read_id = ? ORDER BY run_number DESC LIMIT 1", (row["id"],)).fetchone()
            item.update(_analysis_report(latest) if latest else {"alignment_id":None,"variants":[],"evidence":{}})
            output.append(item)
    return output


@app.get("/api/sanger-reads/{read_id}/analyses")
def list_sanger_analyses(read_id: str, limit: int = Query(default=10, ge=1, le=100), offset: int = Query(default=0, ge=0)) -> dict:
    with connect() as connection:
        connection.execute("BEGIN")
        if not connection.execute("SELECT id FROM sequencing_reads WHERE id = ?", (read_id,)).fetchone():
            raise HTTPException(status_code=404, detail="Sanger read not found")
        total = connection.execute("SELECT COUNT(*) FROM sanger_analysis_runs WHERE sequencing_read_id = ?", (read_id,)).fetchone()[0]
        rows = connection.execute("SELECT * FROM sanger_analysis_runs WHERE sequencing_read_id = ? ORDER BY run_number DESC LIMIT ? OFFSET ?", (read_id, limit, offset)).fetchall()
    return {"items":[_analysis_report(row) for row in rows],"total":total,"offset":offset,"limit":limit,"has_more":offset+len(rows)<total}


@app.get("/api/sanger-reads/{read_id}/trace")
def get_sanger_trace(
    read_id: str, start: int = Query(default=0, ge=0), count: int = Query(default=24, ge=1, le=80),
) -> dict[str, object]:
    with connect() as connection:
        row = connection.execute(
            "SELECT storage_path, file_sha256 FROM sequencing_reads WHERE id = ?", (read_id,)
        ).fetchone()
    if not row:
        raise HTTPException(status_code=404, detail="Sanger read not found")
    path = Path(row["storage_path"])
    try:
        if file_sha256(path) != row["file_sha256"]:
            raise HTTPException(status_code=409, detail="Original AB1 hash changed; trace evidence cannot be trusted")
        result = trace_window(path, start, count)
    except OSError as exc:
        raise HTTPException(status_code=409, detail="Original AB1 file is unavailable") from exc
    except SangerReadError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return {"read_id": read_id, "file_sha256": row["file_sha256"], **result}
