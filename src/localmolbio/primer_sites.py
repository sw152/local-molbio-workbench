"""Exact oligo matches on a single reference, not a specificity prediction."""
from Bio.Seq import Seq


def exact_reference_sites(template: str, oligo: str, topology: str, limit: int = 100) -> dict:
    sequence, primer = template.upper(), oligo.upper()
    if not sequence or not primer or set(sequence + primer) - set("ACGT"):
        raise ValueError("Exact-site review requires nonempty A/C/G/T sequences")
    if len(primer) > len(sequence) or not 1 <= limit <= 1000:
        raise ValueError("Invalid primer length or stored-site limit")
    if topology not in {"linear", "circular", "unknown"}:
        raise ValueError("Unsupported reference topology")
    length = len(sequence)
    haystack = sequence + (sequence[:len(primer)-1] if topology == "circular" else "")
    sites, total = [], 0
    # Palindromic oligos have two directional matches at the same interval.
    for strand, needle in ((1, primer), (-1, str(Seq(primer).reverse_complement()))):
        position = haystack.find(needle)
        while position >= 0:
            total += 1
            if len(sites) < limit:
                stop = position + len(primer)
                wraps = stop > length
                sites.append({"start": position, "end": stop - length if wraps else stop,
                              "strand": strand, "wraps_origin": wraps,
                              "segments": [[position, length], [0, stop-length]] if wraps else [[position, stop]]})
            position = haystack.find(needle, position + 1)
    return {"method": "full_length_exact_match_v1", "scope": "current_reference_only",
            "reference_length": length, "reference_topology": topology,
            "origin_checked": topology == "circular", "mismatches_checked": False,
            "coordinate_system": "zero-based-half-open", "total_directional_matches": total,
            "sites": sites, "stored_site_limit": limit, "truncated": total > len(sites)}
