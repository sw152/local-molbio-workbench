"""Minimal synthetic ABIF binary; not an instrument trace or a peak-shape fixture."""
import struct


def synthetic_ab1(sequence: str, qualities: list[int]) -> bytes:
    if len(sequence) != len(qualities):
        raise ValueError("Expected one quality per base")
    tags = [(b"PBAS", 2, sequence.encode("ascii")), (b"PCON", 2, bytes(qualities)), (b"SMPL", 1, b"synthetic")]
    directory_offset = 30
    payload_offset = directory_offset + 28 * len(tags)
    entries, payload = [], bytearray()
    for name, number, data in tags:
        offset = int.from_bytes(data.ljust(4, b"\0"), "big") if len(data) <= 4 else payload_offset + len(payload)
        entries.append(struct.pack(">4sI2H4I", name, number, 2, 1, len(data), len(data), offset, 0))
        if len(data) > 4:
            payload.extend(data)
    return b"ABIF" + struct.pack(">H4sI2H3I", 101, b"tdir", 1, 1023, 28, len(tags), 28 * len(tags), directory_offset) + b"".join(entries) + payload
