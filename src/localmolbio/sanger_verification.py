from __future__ import annotations

from dataclasses import dataclass
from itertools import islice

from Bio import __version__ as biopython_version
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


class SangerVerificationError(ValueError):
    """Raised when a read cannot produce a meaningful alignment."""


ALIGNMENT_PARAMETERS = {
    "algorithm": "Biopython.PairwiseAligner",
    "algorithm_version": biopython_version,
    "evidence_version": 2,
    "mode": "local",
    "match_score": 2,
    "mismatch_score": -2,
    "open_gap_score": -8,
    "extend_gap_score": -0.5,
    "wildcard": "N",
    "min_aligned_bases": 12,
    "max_dp_cells": 20_000_000,
    "candidate_limit_per_direction": 64,
    "quality_review_threshold": 20,
}


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
    evidence: dict[str, object]


def align_sanger_read(
    reference: str, read: str, direction: str, circular: bool,
    qualities: list[int] | None = None,
) -> SangerAlignment:
    """Return local read evidence, never a whole-plasmid verification verdict.

    Coordinates are zero-based, half-open; insertions use reference boundaries.
    Query positions refer to the oriented read. Original positions are also saved.
    Identity includes gaps and ambiguous columns in its denominator.
    """
    template, raw_query = reference.upper(), read.upper()
    if direction not in {"forward", "reverse", "unknown"}:
        raise SangerVerificationError("Direction must be forward, reverse or unknown")
    if len(template) < 12 or len(raw_query) < 12:
        raise SangerVerificationError("Reference and read must each contain at least 12 bases")
    if set(template + raw_query) - set("ACGTN"):
        raise SangerVerificationError("Alignment supports only A/C/G/T/N base calls")
    if qualities is not None and (
        len(qualities) != len(raw_query)
        or any(not isinstance(q, int) or q < 0 for q in qualities)
    ):
        raise SangerVerificationError("Quality values must be nonnegative integers, one per base")
    target = template * 2 if circular else template
    if len(target) * len(raw_query) > ALIGNMENT_PARAMETERS["max_dp_cells"]:
        raise SangerVerificationError("Read/reference pair exceeds interactive alignment limit")
    aligner = PairwiseAligner(mode="local")
    for key in ("match_score", "mismatch_score", "open_gap_score", "extend_gap_score", "wildcard"):
        setattr(aligner, key, ALIGNMENT_PARAMETERS[key])
    orientations = ("forward", "reverse") if direction == "unknown" else (direction,)
    candidates = []
    flags = []
    for orientation in orientations:
        query = str(Seq(raw_query).reverse_complement()) if orientation == "reverse" else raw_query
        alignments = aligner.align(target, query)
        # Bounded enumeration avoids exponential traceback on repeated sequences.
        for index, alignment in enumerate(islice(alignments, 65)):
            if index == 64:
                flags.append("candidate_search_truncated")
                break
            start, end = map(int, (alignment.coordinates[0, 0], alignment.coordinates[0, -1]))
            if circular and end - start > len(template):
                continue
            signature = (orientation, tuple(
                (int(t) % len(template) if t >= 0 else -1, int(q))
                for t, q in zip(*alignment.indices)
            ))
            candidates.append((alignment.score, orientation, query, alignment, signature))
    if not candidates:
        raise SangerVerificationError("No single-traversal local alignment found")
    best_score = max(item[0] for item in candidates)
    best = [item for item in candidates if abs(item[0] - best_score) < 1e-6]
    _, orientation, query, alignment, _ = best[0]
    if len({item[4] for item in best}) > 1:
        flags.append("ambiguous_alignment")
    oriented_quality = qualities[::-1] if qualities is not None and orientation == "reverse" else qualities
    variants = []
    matched = mismatched = inserted = deleted = ambiguous = 0
    positions, query_positions, covered_positions = [], [], set()
    previous_reference = None
    low_quality = 0
    for target_index, query_index in zip(*alignment.indices):
        t, q = int(target_index), int(query_index)
        position = t % len(template) if t >= 0 else (previous_reference + 1) % len(template) if circular else previous_reference + 1
        if t >= 0:
            positions.append(t)
            previous_reference = t
        if q >= 0:
            query_positions.append(q)
        quality = oriented_quality[q] if q >= 0 and oriented_quality is not None else None
        if quality is not None and quality < 20:
            low_quality += 1
        if t >= 0 and q >= 0:
            covered_positions.add(position)
            if "N" in (target[t], query[q]):
                ambiguous += 1
                kind = "ambiguous"
            elif target[t] == query[q]:
                matched += 1
                continue
            else:
                mismatched += 1
                kind = "substitution"
        elif t < 0:
            inserted += 1
            kind = "insertion"
        else:
            deleted += 1
            kind = "deletion"
        variants.append({
            "kind": kind, "position": position,
            "reference": target[t] if t >= 0 else "", "read": query[q] if q >= 0 else "",
            "read_position": q if q >= 0 else None,
            "original_read_position": (len(query) - 1 - q if orientation == "reverse" else q) if q >= 0 else None,
            "phred": quality,
        })
    aligned = matched + mismatched + ambiguous
    if aligned < 12:
        raise SangerVerificationError("Local match spans fewer than 12 aligned bases")
    start = positions[0] % len(template)
    span = positions[-1] - positions[0] + 1
    wraps = circular and start + span > len(template)
    end = start + span - len(template) if wraps else start + span
    read_fraction = len(query_positions) / len(query)
    if read_fraction < 0.8:
        flags.append("partial_read_alignment")
    if ambiguous:
        flags.append("ambiguous_bases")
    if qualities is None:
        flags.append("quality_unavailable")
    elif low_quality:
        flags.append("low_quality_aligned_bases")
    intervals = []
    for p in sorted(covered_positions):
        if intervals and intervals[-1][1] == p:
            intervals[-1][1] = p + 1
        else:
            intervals.append([p, p + 1])
    evidence = {
        "parameters": ALIGNMENT_PARAMETERS.copy(),
        "coordinate_system": "zero-based-half-open; insertion positions are boundaries",
        "variant_representation": "per-base, not normalized",
        "direction": orientation, "requested_direction": direction, "score": best_score,
        "read_start": min(query_positions), "read_end": max(query_positions) + 1,
        "read_aligned_fraction": round(read_fraction, 6),
        "reference_covered_intervals": intervals,
        "reference_covered_fraction": round(len(covered_positions) / len(template), 6),
        "ambiguous_bases": ambiguous, "low_quality_aligned_bases": low_quality if qualities is not None else None,
        "review_flags": sorted(set(flags)), "scope": "local_read_alignment",
        "whole_reference_verified": False,
    }
    return SangerAlignment(
        start, end, wraps, aligned, matched, mismatched, inserted, deleted,
        round(matched / (aligned + inserted + deleted), 6), variants, evidence,
    )
