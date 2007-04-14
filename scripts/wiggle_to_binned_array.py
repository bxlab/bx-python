#!/usr/bin/env python

"""
Convert wiggle data to a binned array. This assumes the input data is on a
single chromosome and does no sanity checks!

usage: %prog score_file out_file < wiggle_data
    -c, --comp=type: compression type (none, zlib, lzo)
"""

from __future__ import division

import sys
import psyco_full
import bx.wiggle
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.cookbook import doc_optparse
from bx import misc

def main():
    
    # Parse command line
    options, args = doc_optparse.parse( __doc__ )
    try:
        if options.comp:
            comp_type = options.comp
        else:
            comp_type = None
        score_fname = args[0]
        out_fname = args[1]
    except:
        doc_optparse.exit()

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
    if comp_type:
        scores.to_file( out, comp_type=comp_type )
    else:    
        scores.to_file( out )
    out.close()

if __name__ == "__main__": main()
