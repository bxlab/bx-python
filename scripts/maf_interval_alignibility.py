#!/usr/bin/env python2.4

"""
WARNING: bz2/bz2t support and file cache support are new and not as well
         tested. 

usage: %prog maf_files [options] < interval_file
    -s, --species=SPECIES: Comma separated list of species to include
    -p, --prefix=PREFIX: Prefix to add to each interval chrom (usually reference species)
   -C, --usecache:   Use a cache that keeps blocks of the MAF files in memory (requires ~20MB per MAF)
"""

from __future__ import division

import psyco_full

from bx.cookbook import doc_optparse

import bx.align.maf
from bx import misc
import os
import sys

from numpy import *

def main():
    # Parse Command Line
    options, args = doc_optparse.parse( __doc__ )
    try:
        maf_files = args
        species = options.species.split( "," )
        prefix = options.prefix
        use_cache = bool( options.usecache )
        if not prefix:
            prefix = ""
    except:
        doc_optparse.exit()
    # Open indexed access to mafs
    index = bx.align.maf.MultiIndexed( maf_files, 
                                      parse_e_rows=True,
                                      use_cache=use_cache )
    # Print header
    print "#chr", "start", "end",
    for s in species:
        print s,
    print
    # Iterate over input ranges 
    for line in sys.stdin:
        fields = line.split()
        # Input is BED3+
        chr, start, end = fields[0], int( fields[1] ), int( fields[2] )
        length = end - start
        assert length > 0, "Interval has length less than one"
        # Prepend prefix if specified
        src = prefix + chr    
        # Keep a bitset for each species noting covered pieces
        aligned_bits = []
        missing_bits = []
        for s in species:
            aligned_bits.append( zeros( length, dtype=bool ) )
            missing_bits.append( zeros( length, dtype=bool ) )
        # Find overlap with reference component
        blocks = index.get( src, start, end )
        # Determine alignability for each position
        for block in blocks:
            # Determine the piece of the human interval this block covers, 
            # relative to the start of the interval of interest
            ref = block.get_component_by_src( src )
            assert ref.strand == "+", \
                "Reference species blocks must be on '+' strand"
            rel_start = max( start, ref.start ) - start
            rel_end = min( end, ref.end ) - start
            # Check alignability for each species
            for i, s in enumerate( species ):
                other = block.get_component_by_src_start( s )
                # Species does not appear at all indicates unaligned (best we
                # can do here?)
                if other is None:
                    continue
                # An empty component might indicate missing data, all other
                # cases (even contiguous) we count as not aligned
                if other.empty:
                    if other.synteny_empty == bx.align.maf.MAF_MISSING_STATUS:
                        missing_bits[i][rel_start:rel_end] = True
                # Otherwise we have a local alignment with some text, call
                # it aligned
                else:
                    aligned_bits[i][rel_start:rel_end] = True
        # Now determine the total alignment coverage of each interval
        print chr, start, end,
        for i, s in enumerate( species ):
            aligned = sum( aligned_bits[i] )
            missing = sum( missing_bits[i] )
            # An interval will be called missing if it is < 100bp and <50% 
            # present, or more than 100bp and less that 50bp present (yes,
            # arbitrary)
            is_missing = False
            if length < 100 and missing > ( length / 2 ):
                print "NA",
            elif length >= 100 and missing > 50:
                print "NA",
            else:
                print aligned / ( length - missing ),
                
        print
         
    # Close MAF files
    index.close()

if __name__ == "__main__": 
    main()
