#!/usr/bin/env python

"""
usage: %prog score_file 

Writes compressed data from a wiggle file by chromosome.
"""

from __future__ import division

import sys
import psyco_full
import bx.wiggle
from bx.binned_array import BinnedArray
from fpconst import isNaN
import cookbook.doc_optparse
import misc


def main():
    
    # Parse command line
    options, args = cookbook.doc_optparse.parse( __doc__ )
    try:
        score_fname = args[0]
    except:
        cookbook.doc_optparse.exit()

    scores = {}
    for i, ( chrom, pos, val ) in enumerate( bx.wiggle.Reader( open(sys.argv[1]) ) ):
        if not chrom in scores: scores[ chrom ] = BinnedArray()
        scores[chrom][pos] = val

        # Status
        if i % 10000 == 0: print i, "scores processed"

    for chr in scores.keys():
        out = open( chr, "w" )
        scores[chr].to_file( out )
        out.close()

if __name__ == "__main__": main()
