#!/usr/bin/env python

"""
Read two lists of intervals (with chromosomes) and count the number of entries
in the second set that intersect any entry in the first set.

TODO: This could use bitsets rather than the intervals package, would it be
      faster?

usage: %prog bed1 bed2 > out
"""

from __future__ import division

import psyco_full

from bx import intervals
from bx import misc
import string
import sys

def main():

    intersecters = {}

    # Read ranges

    for chr, start, end in read_intervals( misc.open_compressed( sys.argv[1] ) ):
        if not intersecters.has_key( chr ): intersecters[ chr ] = intervals.Intersecter()
        intersecters[ chr ].add_interval( intervals.Interval( start, end ) )

    # Count intersection

    total = 0

    for chr, start, end in read_intervals( misc.open_compressed( sys.argv[2] ) ):
        if intersecters.has_key( chr ):
            intersection = intersecters[ chr ].find( start, end )
            if intersection: 
                #print chr, intersection
                total += 1

    print total

def read_intervals( input ):
    for line in input:
        fields = line.split()
        yield fields[0], int( fields[1] ), int( fields[2] )

main()
