#!/usr/bin/env python
"""
Take an input of GenomicIntervalReaders, subtract, and yield
GenomicIntervals.  The returned GenomicIntervals will be in the order
of the first set of intervals passed in, with the corresponding
meta-data.
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def subtract(readers, mincols=1, upstream_pad=0, downstream_pad=0, pieces=True, lens={}, comments=True):

    # Read all but first into bitsets and union to one (if confused,
    # read DeMorgan's...)
    primary = readers[0]
    union = readers[1:]
    bitsets = union[0].binned_bitsets(upstream_pad = upstream_pad, downstream_pad = downstream_pad, lens = lens)
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
        if type( interval ) is Header:
            yield interval
        if type( interval ) is Comment and comments:
            yield interval
        elif type( interval ) == GenomicInterval:
            chrom = interval.chrom
            start = int(interval.start)
            end = int(interval.end)
            if chrom not in bitsets: continue
            if start > end: warn( "Interval start after end! on line '%d' of second input" % f.lineno() )
            out_intervals = []
            # Find the intervals that meet the criteria (for the three sensible
            # permutations of reverse and pieces)
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
