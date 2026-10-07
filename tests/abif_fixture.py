"""Synthetic ABIF binary and optional artificial peaks; never instrument data."""
import struct
import math


def synthetic_ab1(sequence: str, qualities: list[int], *, trace: bool = False, order: str = "GATC", peak_override=None) -> bytes:
    if len(sequence) != len(qualities):
        raise ValueError("Expected one quality per base")
    tags = [(b"PBAS", 2, 2, 1, sequence.encode("ascii")), (b"PCON", 2, 2, 1, bytes(qualities)), (b"SMPL", 1, 2, 1, b"synthetic")]
    if trace:
        peaks = [12 + i * 12 for i in range(len(sequence))]
        locations = peaks if peak_override is None else peak_override
        tags += [(b"FWO_", 1, 2, 1, order.encode("ascii")), (b"PLOC", 2, 4, 2, struct.pack(">" + "h" * len(locations), *locations))]
        for i, base in enumerate(order):
            values = [int(sum((900 if sequence[j] == base else 40) * math.exp(-((sample-p)/2.2)**2) for j, p in enumerate(peaks) if abs(sample-p) < 10)) for sample in range(peaks[-1]+13)]
            tags.append((b"DATA", 9+i, 4, 2, struct.pack(">" + "h" * len(values), *values)))
    directory_offset = 30
    payload_offset = directory_offset + 28 * len(tags)
    entries, payload = [], bytearray()
    for name, number, code, size, data in tags:
        offset = int.from_bytes(data.ljust(4, b"\0"), "big") if len(data) <= 4 else payload_offset + len(payload)
        entries.append(struct.pack(">4sI2H4I", name, number, code, size, len(data)//size, len(data), offset, 0))
        if len(data) > 4:
            payload.extend(data)
    return b"ABIF" + struct.pack(">H4sI2H3I", 101, b"tdir", 1, 1023, 28, len(tags), 28 * len(tags), directory_offset) + b"".join(entries) + payload
