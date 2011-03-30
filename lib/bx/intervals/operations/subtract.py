#!/usr/bin/env python
"""
Subtract one set of genomic intervals from another (base-by-base or whole
intervals). The returned GenomicIntervals will be in the order
of the first set of intervals passed in, with the corresponding
meta-data.
"""

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def subtract(readers, mincols=1, upstream_pad=0, downstream_pad=0, pieces=True, lens={}, comments=True):
    # The incoming lens dictionary is a dictionary of chromosome lengths which are used to initialize the bitsets.
    # Read all but first into bitsets and union to one (if confused, read DeMorgan's...)
    primary = readers[0]
    union = readers[1:]
    # Handle any ValueError, IndexError and OverflowError exceptions that may be thrown when
    # the bitsets are being created by skipping the problem lines
    union[0] = BitsetSafeReaderWrapper( union[0], lens=lens )
    bitsets = union[0].binned_bitsets( upstream_pad=upstream_pad, downstream_pad=downstream_pad, lens=lens )
    union = union[1:]
    for andset in union:
        bitset2 = andset.binned_bitsets(upstream_pad = upstream_pad, downstream_pad = downstream_pad, lens = lens)
        for chrom in bitset2:
            if chrom not in bitsets:
                bitsets[chrom] = bitset2[chrom]
            else:
                bitsets[chrom].ior(bitset2[chrom])
    
    # Read remaining intervals and subtract
    for interval in primary:
        if isinstance(interval, Header):
            yield interval
        if isinstance(interval, Comment) and comments:
            yield interval
        elif isinstance(interval, GenomicInterval):
            chrom = interval.chrom
            if chrom not in bitsets:
                yield interval
            else:
                start = int(interval.start)
                end = int(interval.end)
                if start > end: warn( "Interval start after end!" )
                out_intervals = []
                # Find the intervals that meet the criteria (for the three sensible
                # permutations of reverse and pieces)
                try:
                    if bitsets[ chrom ].count_range( start, end-start ) >= mincols:                
                        if pieces:
                            out_intervals = bits_clear_in_range( bitsets[chrom], start, end )
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
                        # This will work only if primary is a NiceReaderWrapper
                        primary.skipped += 1
                        # no reason to stuff an entire bad file into memmory
                        if primary.skipped < 10:
                            primary.skipped_lines.append( ( primary.linenum, primary.current_line, str( e ) ) )
                    except:
                        pass
                    continue
