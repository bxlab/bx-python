#!/usr/bin/env python
"""
Take an input of GenomicIntervalReaders, intersect, and return a
GenomicIntervalReader.  The returned GenomicIntervalReader will be in
the order of the first set of intervals passed in, with the
corresponding meta-data.

The parameter output can be any object with an append(foo) method.  A
list is used by default, but GenomicIntervalWriter could be used as
well (not yet written as of 3/26/06).
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def intersect(readers, mincols=1, upstream_pad=0, downstream_pad=0, pieces=True, lens={}, comments=True):

    # Read all but first into bitsets and intersect to one
    primary = readers[0]
    intersect = readers[1:]
    bitsets = intersect[0].binned_bitsets(upstream_pad = upstream_pad, downstream_pad = downstream_pad, lens = lens)
    intersect = intersect[1:]
    for andset in intersect:
        bitset2 = andset.binned_bitsets(upstream_pad = upstream_pad, downstream_pad = downstream_pad, lens = lens)
        for chrom in bitsets:
            if chrom not in bitset2: continue
            bitsets[chrom].iand(bitset2[chrom])
        intersect = intersect[1:]
    
    # Read remaining intervals and intersect
    for interval in primary:
        if type( interval ) is Header:
            yield interval
        if type( interval ) is Comment and comments:
            yield interval
        elif type( interval ) == GenomicInterval:
            chrom = interval["chrom"]
            start = int(interval["start"])
            end = int(interval["end"])
            if chrom not in bitsets: continue
            if start > end: warn( "Interval start after end! on line '%d' of second input" % f.lineno() )
            out_intervals = []
            # Intersect or Overlap
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
