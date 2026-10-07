"""Bounded FASTQ inspection with explicit Phred+33 interpretation; no alignment."""
from collections import Counter
from dataclasses import dataclass
import gzip
from pathlib import Path
import zlib

DNA = frozenset(b'ACGTRYSWKMBDHVN')


class FastqError(ValueError):
    pass


@dataclass(frozen=True)
class FastqLimits:
    raw_bytes: int = 2 * 1024**3
    decoded_bytes: int = 8 * 1024**3
    record_bases: int = 4 * 1024**2
    records: int = 10_000_000

    def __post_init__(self):
        if any(type(v) is not int or v <= 0 for v in vars(self).values()):
            raise ValueError('FASTQ limits must be positive integers')


def inspect_fastq(path: Path, compressed: bool, quality_encoding: str, limits=FastqLimits()):
    if quality_encoding != 'phred33':
        raise FastqError('quality_encoding_must_be_explicit_phred33')
    if path.stat().st_size > limits.raw_bytes:
        raise FastqError('raw_size_limit_exceeded')
    decoded = records = bases = ambiguous = n_bases = 0
    min_length = max_length = None
    quality_counts = Counter()
    try:
        with path.open('rb') as raw:
            is_gzip = raw.read(2) == b'\x1f\x8b'
            raw.seek(0)
            if is_gzip != compressed:
                raise FastqError('compression_does_not_match_filename')
            stream = gzip.GzipFile(fileobj=raw, mode='rb') if compressed else raw
            try:
                def line():
                    nonlocal decoded
                    data = stream.readline(limits.record_bases + 3)
                    decoded += len(data)
                    if decoded > limits.decoded_bytes:
                        raise FastqError('decoded_size_limit_exceeded')
                    if not data:
                        return None
                    if len(data) > limits.record_bases + 2:
                        raise FastqError('line_size_limit_exceeded')
                    if data.endswith(b'\n'):
                        data = data[:-1]
                        if data.endswith(b'\r'):
                            data = data[:-1]
                    return data

                while True:
                    header = line()
                    if header is None:
                        break
                    if not header.startswith(b'@') or not header[1:].strip() or any(c != 9 and not 32 <= c <= 126 for c in header):
                        raise FastqError('invalid_fastq_header')
                    sequence_length = sequence_ambiguous = sequence_n = 0
                    while True:
                        part = line()
                        if part is None:
                            raise FastqError('truncated_sequence')
                        if part.startswith(b'+'):
                            if part[1:] and part[1:] != header[1:]:
                                raise FastqError('repeated_header_mismatch')
                            break
                        upper = part.upper()
                        if not upper or not set(upper).issubset(DNA):
                            raise FastqError('invalid_dna_sequence')
                        sequence_length += len(upper)
                        if sequence_length > limits.record_bases:
                            raise FastqError('record_size_limit_exceeded')
                        sequence_n += upper.count(b'N')
                        sequence_ambiguous += len(upper) - sum(upper.count(bytes([c])) for c in b'ACGT')
                    if not sequence_length:
                        raise FastqError('empty_read')
                    quality_length = 0
                    while quality_length < sequence_length:
                        quality = line()
                        if quality is None:
                            raise FastqError('truncated_quality')
                        if not quality or any(c < 33 or c > 126 for c in quality):
                            raise FastqError('invalid_phred33_character')
                        quality_length += len(quality)
                        if quality_length > sequence_length:
                            raise FastqError('sequence_quality_length_mismatch')
                        quality_counts.update(quality)
                    records += 1
                    if records > limits.records:
                        raise FastqError('record_count_limit_exceeded')
                    bases += sequence_length
                    ambiguous += sequence_ambiguous
                    n_bases += sequence_n
                    min_length = sequence_length if min_length is None else min(min_length, sequence_length)
                    max_length = sequence_length if max_length is None else max(max_length, sequence_length)
            finally:
                if compressed:
                    stream.close()
    except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise FastqError('invalid_or_truncated_gzip') from exc
    if not records:
        raise FastqError('empty_fastq')
    return {'schema_version': 1, 'parser': 'bounded-fastq-v1', 'scope': 'fastq_format_and_quality_summary',
            'quality_encoding': 'phred33', 'quality_encoding_source': 'explicit_user_declaration',
            'compression': 'gzip' if compressed else 'none', 'records': records, 'bases': bases,
            'min_read_length': min_length, 'max_read_length': max_length, 'mean_read_length': bases / records,
            'n_bases': n_bases, 'ambiguous_bases': ambiguous,
            'min_phred': min(quality_counts) - 33, 'max_phred': max(quality_counts) - 33,
            'mean_phred': sum((q - 33) * n for q, n in quality_counts.items()) / bases,
            'q20_bases': sum(n for q, n in quality_counts.items() if q >= 53),
            'q30_bases': sum(n for q, n in quality_counts.items() if q >= 63),
            'decoded_bytes': decoded, 'limits': vars(limits), 'read_layout': 'single_file_unpaired',
            'analysis_performed': False, 'whole_reference_verified': False}
