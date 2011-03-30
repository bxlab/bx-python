"""
Compute the intersection of two sets of genomic intervals, either base-by-base
or at the interval level. The returned GenomicIntervalReader will be in
the order of the first set of intervals passed in, with the corresponding 
additional fields.
"""

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def intersect(readers, mincols=1, upstream_pad=0, downstream_pad=0, pieces=True, lens={}, comments=True):
    # The incoming lens dictionary is a dictionary of chromosome lengths which are used to initialize the bitsets.
    # Read all but first into bitsets and intersect to one
    primary = readers[0]
    intersect = readers[1:]
    # Handle any ValueError, IndexError and OverflowError exceptions that may be thrown when
    # the bitsets are being created by skipping the problem lines
    intersect[0] = BitsetSafeReaderWrapper( intersect[0], lens=lens )
    bitsets = intersect[0].binned_bitsets( upstream_pad=upstream_pad, downstream_pad=downstream_pad, lens=lens )
    intersect = intersect[1:]
    for andset in intersect:
        bitset2 = andset.binned_bitsets(upstream_pad = upstream_pad, downstream_pad = downstream_pad, lens = lens)
        for chrom in bitsets:
            if chrom not in bitset2: continue
            bitsets[chrom].iand(bitset2[chrom])
        intersect = intersect[1:]
    
    # Read remaining intervals and intersect
    for interval in primary:
        if isinstance(interval, Header):
            yield interval
        if isinstance(interval, Comment) and comments:
            yield interval
        elif isinstance(interval, GenomicInterval):
            chrom = interval.chrom
            start = int( interval.start )
            end = int( interval.end )
            if chrom not in bitsets:
                continue
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
            out_intervals = []
            # Intersect or Overlap
            try:
                if bitsets[ chrom ].count_range( start, end-start ) >= mincols:                
                    if pieces:
                        out_intervals = bits_set_in_range( bitsets[chrom], start, end )
                    else:
                        out_intervals = [ ( start, end ) ]
                # Write the intervals
                for start, end in out_intervals:
                    new_interval = interval.copy()
                    new_interval.start = start
                    new_interval.end = end
                    yield new_interval
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
