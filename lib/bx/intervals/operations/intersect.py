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

def intersect(intervals, mincols=1, upstream_pad=0, downstream_pad=0, pieces=True, output=[], lens={}, comments=True):

    # Read all but first into bitsets and intersect to one
    primary = intervals[0]
    intersect = intervals[1:]
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
        if type( interval ) is Header: continue
        if type( interval ) is Comment and comments:
            output.append(interval)
        elif type( interval ) == GenomicInterval:
            chrom = interval["chrom"]
            start = int(interval["start"])
            end = int(interval["end"])
            if chrom not in bitsets: continue
            if start > end: warn( "Interval start after end! on line '%d' of second input" % f.lineno() )
            out_intervals = []
            # Find the intervals that meet the criteria (for the three sensible
            # permutations of reverse and pieces)
            if bitsets[ chrom ].count_range( start, end-start ) >= mincols:                
                if pieces:
                    out_intervals = bits_set_in_range( bitsets[chrom], start, end )
                else:
                    out_intervals = [ ( start, end ) ]
            # Write the intervals
            for start, end in out_intervals:
                new_interval = list( interval )
                new_interval[interval.start_col] = str( start )
                new_interval[interval.end_col] = str( end )
                output.append(new_interval)
    return output

# test function
def main():
    f = fileinput.FileInput( "test1.bed" )
    g1 = GenomicIntervalReader( f )
    f2 = fileinput.FileInput( "test2.bed" )
    g2 = GenomicIntervalReader( f2 )
    f3 = fileinput.FileInput( "test3.bed" )
    g3 = GenomicIntervalReader( f3 )
    output = intersect([g1, g2, g3])
    for interval in output:
        print interval

main()
