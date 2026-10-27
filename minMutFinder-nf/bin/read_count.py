#!/usr/bin/env python3

import gzip


def read_count(R1f, R2f):
    """Count total reads in paired FASTQ files.

    Accepts either file paths (str) or already-open file-like objects.
    Handles both plain and gzip-compressed FASTQ.
    """
    def _count(stream):
        total = 0
        for _ in stream:
            total += 1
        return total // 4

    def _iter_lines(path_or_stream):
        if isinstance(path_or_stream, str):
            if path_or_stream.endswith('.gz'):
                return gzip.open(path_or_stream, 'rt')
            return open(path_or_stream, 'r')
        return path_or_stream

    r1_stream = _iter_lines(R1f)
    r2_stream = _iter_lines(R2f)

    r1r = _count(r1_stream)
    r2r = _count(r2_stream)

    rt = r1r + r2r
    return int(rt)