"""Windowed, unmodified analyzed ABIF traces in original uploaded-read order.

Use the unedited PBAS2/PCON2/PLOC2 set together and map DATA9..12 with FWO_1.
Never infer a missing dye order, substitute raw channels, or synthesize traces.
"""
from pathlib import Path

from Bio import SeqIO

from .sanger import SangerReadError


def trace_window(path: Path, start: int = 0, count: int = 24) -> dict:
    if start < 0 or not 1 <= count <= 80:
        raise SangerReadError("Trace window must start at a nonnegative base and contain 1–80 bases")
    try:
        record = SeqIO.read(str(path), "abi")
    except Exception as exc:
        raise SangerReadError("Could not read the original AB1 chromatogram") from exc
    sequence = str(record.seq).upper()
    if not sequence or start >= len(sequence):
        raise SangerReadError("Trace window starts outside the read")
    raw = record.annotations.get("abif_raw", {})
    required = ("FWO_1", "PLOC2", "PCON2", "PBAS2", "DATA9", "DATA10", "DATA11", "DATA12")
    missing = [tag for tag in required if tag not in raw]
    if missing:
        return {"available": False, "reason": "Analyzed chromatogram tags are missing: " + ", ".join(missing)}
    order = raw["FWO_1"]
    if isinstance(order, bytes):
        order = order.decode("ascii", errors="replace")
    if not isinstance(order, str) or len(order) != 4 or set(order) != set("ACGT"):
        raise SangerReadError("AB1 dye order must identify exactly A, C, G and T")
    qualities = record.letter_annotations.get("phred_quality", [])
    peaks = raw["PLOC2"]
    channels = {base: raw[f"DATA{9 + i}"] for i, base in enumerate(order)}
    if not isinstance(peaks, (list, tuple)) or len(peaks) != len(sequence) or len(qualities) != len(sequence):
        raise SangerReadError("AB1 base, peak-location and quality counts do not match")
    if any(not isinstance(channel, (list, tuple)) for channel in channels.values()):
        raise SangerReadError("AB1 analyzed channels are not sample arrays")
    sizes = {len(channel) for channel in channels.values()}
    if len(sizes) != 1 or not next(iter(sizes)):
        raise SangerReadError("AB1 analyzed channels have inconsistent lengths")
    sample_count = next(iter(sizes))
    if any(not isinstance(p, int) or p < 0 or p >= sample_count for p in peaks) or any(a >= b for a, b in zip(peaks, peaks[1:])):
        raise SangerReadError("AB1 peak locations must increase within the analyzed channels")
    end = min(len(sequence), start + count)
    left = (peaks[start - 1] + peaks[start]) // 2 if start else 0
    right = (peaks[end - 1] + peaks[end]) // 2 + 1 if end < len(peaks) else sample_count
    if right - left > 20_000:
        raise SangerReadError("Trace window exceeds 20,000 samples; request fewer bases")
    window = {base: list(channel[left:right]) for base, channel in channels.items()}
    if any(not isinstance(value, int) for channel in window.values() for value in channel):
        raise SangerReadError("AB1 analyzed channels contain invalid samples")
    return {
        "available": True, "orientation": "original_uploaded_read", "read_length": len(sequence),
        "base_start": start, "base_end": end, "sample_start": left, "sample_end": right,
        "source_channel_order": order, "channels": window,
        "bases": [{"position": i, "base": sequence[i], "phred": int(qualities[i]), "sample": peaks[i]} for i in range(start, end)],
        "processing": "stored analyzed signal; no smoothing, resampling or individual-channel normalization",
    }
