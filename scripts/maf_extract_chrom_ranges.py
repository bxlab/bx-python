#!/usr/bin/env python

"""
Reads a list of intervals and a maf. Produces a new maf containing the
portions of the original that overlapped the intervals

NOTE: See maf_extract_ranges_indexed.py which works better / faster for many
      use cases.

TODO: Combine with maf_extract_ranges, and possibly share some code with 
      maf_extract_ranges_indexed.

usage: %prog interval_file refname|refindex [options] < maf_file
   -m, --mincols=10: Minimum length (columns) required for alignment to be output
   -p, --prefix=PREFIX: Prefix
"""

import psyco_full

from bx.cookbook import doc_optparse

import bx.align.maf
from bx import intervals
import sys


def __main__():

    # Parse Command Line

    options, args = doc_optparse.parse( __doc__ )

    try:
        range_filename = args[ 0 ]
        try: 
            refindex = int( args[ 1 ] )
            refname = None
        except: 
            refindex = None
            refname = args[ 1 ]
        if options.mincols: mincols = int( options.mincols )
        else: mincols = 10
        if options.prefix: prefix = options.prefix
        else: prefix = ""
    except:
        doc_optparse.exit()

    # Load Intervals

    intersecters = dict()    
    for line in file( range_filename ):
        fields = line.split()
        src = prefix + fields[0]
        if not src in intersecters: intersecters[src] = intervals.Intersecter()
        intersecters[src].add_interval( intervals.Interval( int( fields[1] ), int( fields[2] ) ) )

    # Start MAF on stdout

    out = bx.align.maf.Writer( sys.stdout )

    # Iterate over input MAF

    for maf in bx.align.maf.Reader( sys.stdin ):
        if refname: 
            sourcenames = [ cmp.src.split('.')[0] for cmp in maf.components ]
            try: refindex = sourcenames.index( refname )
            except:
                continue

        ref_component = maf.components[ refindex ]
        # Find overlap with reference component
        if not ( ref_component.src in intersecters ): continue
        intersections = intersecters[ ref_component.src ].find( ref_component.start, ref_component.end )
        # Keep output maf ordered
        intersections.sort()
        # Write each intersecting block
        for interval in intersections: 
            start = max( interval.start, ref_component.start )
            end = min( interval.end, ref_component.end )
            sliced = maf.slice_by_component( refindex, start, end ) 
            good = True
            for c in sliced.components: 
                if c.size < 1: 
                    good = False
            if good and sliced.text_size > mincols: out.write( sliced )
         
    # Close output MAF

    out.close()

if __name__ == "__main__": __main__()
