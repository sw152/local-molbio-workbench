from __future__ import annotations

from contextlib import asynccontextmanager
import json
from pathlib import Path
from tempfile import NamedTemporaryFile

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse

from .config import database_path
from .database import connect, initialise
from .importer import ImportErrorDetail, _parse_genbank, import_benchling_archive


@asynccontextmanager
async def lifespan(_: FastAPI):
    initialise()
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


@app.get("/api/sequences/{sequence_id}/map")
def sequence_map(sequence_id: str) -> dict[str, object]:
    with connect() as connection:
        row = connection.execute(
            """
            SELECT id, display_name, length_bp, topology, feature_count,
                   features_json, parse_warning_count, parse_warnings_json,
                   archive_member_name, raw_genbank
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
