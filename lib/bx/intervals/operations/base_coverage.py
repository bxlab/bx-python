"""
Determine the number of bases covered by a set of intervals.
"""

from bx.intervals.io import BitsetSafeReaderWrapper
from bx.intervals.operations import MAX_END


def base_coverage(reader):
    # Handle any ValueError, IndexError and OverflowError exceptions that may be thrown when
    # the bitsets are being created by skipping the problem lines
    base_reader = BitsetSafeReaderWrapper(reader, lens={})
    bitsets = base_reader.binned_bitsets()
    coverage = 0
    for chrom in bitsets:
        try:
            coverage += bitsets[chrom].count_range(0, MAX_END)
        except IndexError as e:
            base_reader.skipped += 1
            # no reason to stuff an entire bad file into memmory
            if base_reader.skipped < 10:
                base_reader.skipped_lines.append((base_reader.linenum, base_reader.current_line, str(e)))
            continue
    return coverage
