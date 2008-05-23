"""
Merge overlapping regions in two sets of genomic intervals.
"""

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

# sorting could make this a less memory intensive operation(?)
def merge( interval, mincols=1 ):
    # Handle any ValueError, IndexError and OverflowError exceptions that may be thrown when
    # the bitsets are being created by skipping the problem lines
    interval = BitsetSafeReaderWrapper( interval, lens={} )
    bitsets = interval.binned_bitsets()
    if interval.header:
        yield interval.header
    for chrom in bitsets:
        bitset = bitsets[chrom]
        output = ["."] * (max(interval.chrom_col, interval.start_col, interval.end_col) + 1)
        output[interval.chrom_col] = chrom
        try:
            for start, end in bits_set_in_range(bitset,0, MAX_END):
                output[interval.start_col] = str(start)
                output[interval.end_col] = str(end)
                yield output
        except IndexError, e:
            try:
                # This will work only if interval is a NiceReaderWrapper
                interval.skipped += 1
                # no reason to stuff an entire bad file into memmory
                if interval.skipped < 10:
                    interval.skipped_lines.append( ( interval.linenum, interval.current_line, str( e ) ) )
            except:
                pass
            continue
