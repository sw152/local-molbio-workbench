from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


class SangerVerificationError(ValueError):
    """Raised when a read cannot produce a meaningful alignment."""


@dataclass(frozen=True)
class SangerAlignment:
    reference_start: int
    reference_end: int
    wraps_origin: bool
    aligned_bases: int
    matched_bases: int
    mismatched_bases: int
    inserted_bases: int
    deleted_bases: int
    identity_fraction: float
    variants: list[dict[str, object]]


def align_sanger_read(reference: str, read: str, direction: str, circular: bool) -> SangerAlignment:
    template = reference.upper()
    query = read.upper()
    if direction == "reverse":
        query = str(Seq(query).reverse_complement())
    if len(template) < 12 or len(query) < 12:
        raise SangerVerificationError("Reference and read must each contain at least 12 bases")
    if set(template + query) - {"A", "C", "G", "T", "N"}:
        raise SangerVerificationError("Alignment supports only IUPAC A/C/G/T/N base calls")

    target = template * 2 if circular else template
    aligner = PairwiseAligner(mode="local")
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5
    alignment = next(aligner.align(target, query), None)
    if alignment is None or alignment.score <= 0:
        raise SangerVerificationError("No meaningful local alignment found")

    reference_indices, read_indices = alignment.indices
    variants: list[dict[str, object]] = []
    matched = mismatched = inserted = deleted = 0
    reference_positions = []
    for reference_index, read_index in zip(reference_indices, read_indices):
        if reference_index >= 0:
            reference_positions.append(int(reference_index) % len(template))
        position = int(reference_index) % len(template) if reference_index >= 0 else None
        if reference_index >= 0 and read_index >= 0:
            reference_base = target[int(reference_index)]
            read_base = query[int(read_index)]
            if reference_base == read_base:
                matched += 1
            else:
                mismatched += 1
                variants.append({"kind": "substitution", "position": position, "reference": reference_base, "read": read_base})
        elif reference_index < 0:
            inserted += 1
            variants.append({"kind": "insertion", "position": position, "read": query[int(read_index)]})
        else:
            deleted += 1
            variants.append({"kind": "deletion", "position": position, "reference": target[int(reference_index)]})
    if not reference_positions:
        raise SangerVerificationError("Alignment contains no reference coverage")
    start = reference_positions[0]
    end = (reference_positions[-1] + 1) % len(template)
    wraps = circular and any(a > b for a, b in zip(reference_positions, reference_positions[1:]))
    aligned = matched + mismatched
    return SangerAlignment(start, end, wraps, aligned, matched, mismatched, inserted, deleted, round(matched / aligned, 6) if aligned else 0.0, variants)
