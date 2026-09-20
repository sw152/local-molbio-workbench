from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
from shutil import copyfile
from uuid import uuid4
import warnings
from zipfile import BadZipFile, ZipFile, ZipInfo

from Bio import SeqIO

from .config import archive_dir
from .database import connect, initialise


SUPPORTED_SUFFIXES = {".gb", ".gbk", ".genbank"}


@dataclass(frozen=True)
class ImportSummary:
    import_id: str
    total_members: int
    imported_records: int
    skipped_members: int
    parse_warning_records: int
    parse_warning_count: int
    archive_sha256: str


class ImportErrorDetail(ValueError):
    """Raised when an archive cannot be safely imported."""


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _feature_label(feature: object, index: int) -> str:
    qualifiers = getattr(feature, "qualifiers", {})
    for key in ("label", "gene", "name", "locus_tag", "note"):
        values = qualifiers.get(key)
        if values:
            return str(values[0])
    return f"{getattr(feature, 'type', 'feature')} {index + 1}"


def _serialise_features(record: object) -> list[dict[str, object]]:
    """Return display metadata while keeping raw GenBank as the source of truth."""
    serialised: list[dict[str, object]] = []
    for index, feature in enumerate(getattr(record, "features", [])):
        location = feature.location
        parts = getattr(location, "parts", []) if location is not None else []
        if location is not None and not parts:
            parts = [location]
        segments = [
            {"start": int(part.start), "end": int(part.end), "strand": part.strand}
            for part in parts
        ]
        serialised.append(
            {
                "index": index,
                "type": feature.type,
                "label": _feature_label(feature, index),
                "location_text": str(location) if location is not None else None,
                "segments": segments,
                "renderable": bool(segments),
            }
        )
    return serialised


def _parse_genbank(raw: bytes | str, member_name: str) -> dict[str, object]:
    if isinstance(raw, str):
        text = raw
    else:
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            text = raw.decode("latin-1")

    try:
        with warnings.catch_warnings(record=True) as captured_warnings:
            warnings.simplefilter("always")
            record = SeqIO.read(StringIO(text), "genbank")
    except Exception as exc:  # Biopython exposes several parser exceptions.
        raise ImportErrorDetail(f"{member_name}: invalid GenBank ({exc})") from exc

    sequence = str(record.seq).upper()
    if not sequence:
        raise ImportErrorDetail(f"{member_name}: GenBank record contains no bases")

    topology = str(record.annotations.get("topology", "unknown")).lower()
    if topology not in {"circular", "linear"}:
        topology = "unknown"

    return {
        "display_name": record.name if record.name not in {".", "<unknown name>"} else Path(member_name).stem,
        "length_bp": len(sequence),
        "topology": topology,
        "molecule_type": record.annotations.get("molecule_type"),
        "feature_count": len(record.features),
        "features": _serialise_features(record),
        "sequence": sequence,
        "parse_warnings": [str(item.message) for item in captured_warnings],
        "sequence_sha256": sha256(sequence.encode("ascii")).hexdigest(),
        "raw_genbank": text,
    }


def ensure_initial_revisions(database: Path | None = None) -> int:
    """Backfill one immutable import revision for records created by earlier releases."""
    initialise(database)
    inserted = 0
    with connect(database) as connection:
        rows = connection.execute(
            """
            SELECT sequences.id, sequences.display_name, sequences.sequence_sha256,
                   sequences.topology, sequences.features_json, sequences.raw_genbank,
                   sequences.archive_member_name, sequences.created_at
            FROM sequences
            LEFT JOIN sequence_revisions
              ON sequence_revisions.sequence_id = sequences.id
            WHERE sequence_revisions.id IS NULL
            """
        ).fetchall()
        for row in rows:
            parsed = _parse_genbank(row["raw_genbank"], row["archive_member_name"])
            features_json = row["features_json"]
            if features_json == "[]":
                features_json = json.dumps(parsed["features"])
            connection.execute(
                """
                INSERT INTO sequence_revisions (
                    id, sequence_id, parent_revision_id, revision_number, label,
                    sequence_text, sequence_sha256, topology, features_json,
                    source_kind, created_at
                ) VALUES (?, ?, NULL, 1, ?, ?, ?, ?, ?, 'import', ?)
                """,
                (
                    str(uuid4()),
                    row["id"],
                    row["display_name"],
                    parsed["sequence"],
                    row["sequence_sha256"],
                    row["topology"],
                    features_json,
                    row["created_at"],
                ),
            )
            inserted += 1
    return inserted


def import_benchling_archive(archive: Path, database: Path | None = None) -> ImportSummary:
    """Import a Benchling GenBank ZIP without flattening duplicate member names."""
    archive = archive.resolve()
    if not archive.is_file():
        raise ImportErrorDetail(f"Archive not found: {archive}")

    try:
        with ZipFile(archive) as zipper:
            members = [member for member in zipper.infolist() if not member.is_dir()]
            genbank_members = [
                member
                for member in members
                if Path(member.filename).suffix.lower() in SUPPORTED_SUFFIXES
            ]
            if not genbank_members:
                raise ImportErrorDetail("Archive contains no supported GenBank files")

            archive_hash = _file_sha256(archive)
            initialise(database)
            ensure_initial_revisions(database)
            with connect(database) as connection:
                prior = connection.execute(
                    "SELECT id FROM imports WHERE archive_sha256 = ?", (archive_hash,)
                ).fetchone()
                if prior:
                    raise ImportErrorDetail("This exact archive has already been imported")

            occurrences: Counter[str] = Counter()
            rows: list[tuple[ZipInfo, int, bytes, dict[str, object]]] = []
            errors: list[str] = []
            for index, member in enumerate(members):
                if Path(member.filename).suffix.lower() not in SUPPORTED_SUFFIXES:
                    continue
                raw = zipper.read(member)
                occurrences[member.filename] += 1
                try:
                    parsed = _parse_genbank(raw, member.filename)
                except ImportErrorDetail as exc:
                    errors.append(str(exc))
                    continue
                rows.append((member, index, raw, parsed))
    except BadZipFile as exc:
        raise ImportErrorDetail("File is not a valid ZIP archive") from exc

    if not rows:
        detail = "; ".join(errors[:3]) or "No GenBank records could be parsed"
        raise ImportErrorDetail(detail)

    destination_dir = archive_dir()
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / f"{archive_hash}.zip"
    if not destination.exists():
        copyfile(archive, destination)

    import_id = str(uuid4())
    imported_at = _utc_now()
    parse_warning_records = sum(bool(parsed["parse_warnings"]) for _, _, _, parsed in rows)
    parse_warning_count = sum(len(parsed["parse_warnings"]) for _, _, _, parsed in rows)
    with connect(database) as connection:
        connection.execute(
            """
            INSERT INTO imports (
                id, original_filename, archive_sha256, imported_at, total_members,
                imported_records, skipped_members, parse_warning_records,
                parse_warning_count, archive_path
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                import_id,
                archive.name,
                archive_hash,
                imported_at,
                len(members),
                len(rows),
                len(errors),
                parse_warning_records,
                parse_warning_count,
                str(destination),
            ),
        )
        seen: Counter[str] = Counter()
        for member, index, raw, parsed in rows:
            seen[member.filename] += 1
            sequence_id = str(uuid4())
            connection.execute(
                """
                INSERT INTO sequences (
                    id, import_id, archive_member_index, archive_member_name,
                    archive_member_occurrence, member_sha256, sequence_sha256,
                    display_name, length_bp, topology, molecule_type, feature_count,
                    features_json, parse_warning_count, parse_warnings_json, raw_genbank,
                    created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    sequence_id,
                    import_id,
                    index,
                    member.filename,
                    seen[member.filename],
                    sha256(raw).hexdigest(),
                    parsed["sequence_sha256"],
                    parsed["display_name"],
                    parsed["length_bp"],
                    parsed["topology"],
                    parsed["molecule_type"],
                    parsed["feature_count"],
                    json.dumps(parsed["features"]),
                    len(parsed["parse_warnings"]),
                    json.dumps(parsed["parse_warnings"]),
                    parsed["raw_genbank"],
                    imported_at,
                ),
            )
            connection.execute(
                """
                INSERT INTO sequence_revisions (
                    id, sequence_id, parent_revision_id, revision_number, label,
                    sequence_text, sequence_sha256, topology, features_json,
                    source_kind, created_at
                ) VALUES (?, ?, NULL, 1, ?, ?, ?, ?, ?, 'import', ?)
                """,
                (
                    str(uuid4()),
                    sequence_id,
                    parsed["display_name"],
                    parsed["sequence"],
                    parsed["sequence_sha256"],
                    parsed["topology"],
                    json.dumps(parsed["features"]),
                    imported_at,
                ),
            )

    return ImportSummary(
        import_id=import_id,
        total_members=len(members),
        imported_records=len(rows),
        skipped_members=len(errors),
        parse_warning_records=parse_warning_records,
        parse_warning_count=parse_warning_count,
        archive_sha256=archive_hash,
    )
