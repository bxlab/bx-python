"""
Determine amount of each interval in one set covered by the intervals of 
another set. Adds two columns to the first input, giving number of bases 
covered and percent coverage on the second input.
"""

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def coverage(readers, comments=True):
    # The incoming lens dictionary is a dictionary of chromosome lengths which are used to initialize the bitsets.
    primary = readers[0]
    intersect = readers[1:]
    # Handle any ValueError, IndexError and OverflowError exceptions that may be thrown when
    # the bitsets are being created by skipping the problem lines
    intersect[0] = BitsetSafeReaderWrapper( intersect[0], lens={} )
    bitsets = intersect[0].binned_bitsets()
    intersect = intersect[1:]
    for andset in intersect:
        bitset2 = andset.binned_bitsets()
        for chrom in bitsets:
            if chrom not in bitset2: continue
            bitsets[chrom].ior(bitset2[chrom])
        intersect = intersect[1:]

    # Read remaining intervals and give coverage
    for interval in primary:
        if isinstance(interval, Header):
            yield interval
        if isinstance(interval, Comment) and comments:
            yield interval
        elif isinstance(interval, GenomicInterval):
            chrom = interval.chrom
            start = int(interval.start)
            end = int(interval.end)
            if start > end:
                try:
                    # This will only work if primary is a NiceReaderWrapper
                    primary.skipped += 1
                    # no reason to stuff an entire bad file into memmory
                    if primary.skipped < 10:
                        primary.skipped_lines.append( ( primary.linenum, primary.current_line, "Interval start after end!" ) )
                except:
                    pass
                continue
            if chrom not in bitsets:
                bases_covered = 0
                percent = 0.0
            else:
                try:
                    bases_covered = bitsets[ chrom ].count_range( start, end-start )
                except IndexError, e:
                    try:
                        # This will only work if primary is a NiceReaderWrapper
                        primary.skipped += 1
                        # no reason to stuff an entire bad file into memmory
                        if primary.skipped < 10:
                            primary.skipped_lines.append( ( primary.linenum, primary.current_line, str( e ) ) )
                    except:
                        pass
                    continue
                if (end - start) == 0:
                    percent = 0
                else:
                    percent = float(bases_covered) / float(end - start)
            interval.fields.append(str(bases_covered))
            interval.fields.append(str(percent))
            yield interval
