from __future__ import annotations

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pytest

from localmolbio.database import connect
from localmolbio.importer import ImportErrorDetail, import_benchling_archive


GENBANK = """LOCUS       {name:<16} {length:>7} bp    ds-DNA     circular     01-JAN-2026
DEFINITION  synthetic test record.
FEATURES             Location/Qualifiers
     source          1..{length}
                     /organism=\"synthetic construct\"
ORIGIN
        1 {bases}
//
"""


def make_record(name: str, bases: str) -> str:
    return GENBANK.format(name=name, length=len(bases), bases=bases.lower())


def test_import_preserves_duplicate_archive_member_names(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    archive = tmp_path / "benchling.zip"
    first = make_record("first", "ATGCATGC")
    second = make_record("second", "ATGCAAAC")
    with ZipFile(archive, "w", ZIP_DEFLATED) as zipper:
        zipper.writestr("same-name.gb", first)
        zipper.writestr("same-name.gb", second)
        zipper.writestr("notes.gp", "not a sequence")

    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    database = tmp_path / "runtime" / "workbench.sqlite3"
    result = import_benchling_archive(archive, database)

    assert result.total_members == 3
    assert result.imported_records == 2
    assert result.parse_warning_count == 0
    with connect(database) as connection:
        rows = connection.execute(
            "SELECT display_name, archive_member_occurrence FROM sequences ORDER BY archive_member_index"
        ).fetchall()
    assert [(row["display_name"], row["archive_member_occurrence"]) for row in rows] == [
        ("first", 1),
        ("second", 2),
    ]


def test_same_archive_cannot_be_imported_twice(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    archive = tmp_path / "benchling.zip"
    with ZipFile(archive, "w", ZIP_DEFLATED) as zipper:
        zipper.writestr("sequence.gb", make_record("one", "ATGCATGC"))

    monkeypatch.setenv("MOLBIO_DATA_DIR", str(tmp_path / "runtime"))
    database = tmp_path / "runtime" / "workbench.sqlite3"
    import_benchling_archive(archive, database)
    with pytest.raises(ImportErrorDetail, match="already been imported"):
        import_benchling_archive(archive, database)


def test_import_records_parser_warnings_for_annotation_review(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = tmp_path / "wrapped-feature.zip"
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
    database = tmp_path / "runtime" / "workbench.sqlite3"
    result = import_benchling_archive(archive, database)

    assert result.imported_records == 1
    assert result.parse_warning_records == 1
    assert result.parse_warning_count == 1
    with connect(database) as connection:
        row = connection.execute(
            "SELECT parse_warning_count, parse_warnings_json FROM sequences"
        ).fetchone()
    assert row["parse_warning_count"] == 1
    assert "spans the origin" in row["parse_warnings_json"]
