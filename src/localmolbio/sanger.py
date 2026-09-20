from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from statistics import mean

from Bio import SeqIO


class SangerReadError(ValueError):
    """Raised when an AB1 file cannot produce a usable Sanger read."""


@dataclass(frozen=True)
class SangerRead:
    sequence: str
    quality_summary: dict[str, float | int]
    parser_metadata: dict[str, str]


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_ab1(path: Path) -> SangerRead:
    try:
        record = SeqIO.read(str(path), "abi")
    except Exception as exc:
        raise SangerReadError(f"Could not parse AB1 file: {exc}") from exc

    sequence = str(record.seq).upper()
    if not sequence:
        raise SangerReadError("AB1 file contains no base calls")
    qualities = [int(value) for value in record.letter_annotations.get("phred_quality", [])]
    if len(qualities) != len(sequence):
        raise SangerReadError("AB1 base-call and quality lengths do not match")

    summary: dict[str, float | int] = {
        "mean_phred": round(mean(qualities), 2),
        "min_phred": min(qualities),
        "q20_bases": sum(value >= 20 for value in qualities),
        "q30_bases": sum(value >= 30 for value in qualities),
        "q20_fraction": round(sum(value >= 20 for value in qualities) / len(qualities), 4),
    }
    raw = record.annotations.get("abif_raw", {})
    metadata = {
        key: str(raw[key])
        for key in ("SMPL1", "CMNT1", "PDMF1", "RUND1")
        if key in raw
    }
    return SangerRead(sequence=sequence, quality_summary=summary, parser_metadata=metadata)
