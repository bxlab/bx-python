#!/usr/bin/env python

"""
usage: %prog score_file out_file
"""

from __future__ import division

import sys
import psyco_full
import bx.wiggle
from bx.binned_array import BinnedArray
from fpconst import isNaN
import cookbook.doc_optparse
import misc

def read_scores( f ):
    scores_by_chrom = dict()

    return scores_by_chrom

def main():
    
    # Parse command line
    options, args = cookbook.doc_optparse.parse( __doc__ )
    try:
        score_fname = args[0]
        out_fname = args[1]
    except:
        cookbook.doc_optparse.exit()

    scores = BinnedArray()
    ## last_chrom = None
    for i, ( chrom, pos, val ) in enumerate( bx.wiggle.Reader( misc.open_compressed( score_fname ) ) ):
        #if last_chrom is None: 
        #    last_chrom = chrom
        #else: 
        #    assert chrom == last_chrom, "This script expects a 'wiggle' input on only one chromosome"
        scores[pos] = val
        # Status
        if i % 10000 == 0: print i, "scores processed"

    out = open( out_fname, "w" )
    scores.to_file( out )
    out.close()

if __name__ == "__main__": main()
