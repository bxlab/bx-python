#!/usr/bin/env python2.3

"""

Reads a maf file from stdin and applies the mapping file specified by
`mapping_file` to produce a sequence of integers. Then for each possible word
of length `motif_len` in this integer alphabet print the number of times
that word occurs in the block. 

usage: %prog motif_len mapping_file < maf_file > counts
"""

from __future__ import division

import psyco_full

import bx.align.maf
from bx import seqmapping
import string
import sys
from Numeric import *

def main():

    word_length = int( sys.argv[1] )
    align_count, alpha_map = seqmapping.alignment_mapping_from_file( file( sys.argv[2] ) )

    for maf in bx.align.maf.Reader( sys.stdin ):
        assert len( maf.components ) == align_count
        # Translate alignment to ints
        ints = seqmapping.DNA.translate_list( [ c.text for c in maf.components ] )
        # Apply mapping 
        ints = alpha_map.translate( ints )
        # Count words
        radix = alpha_map.get_out_size()
        counts = zeros( radix ** word_length, Int ) 
        total = 0
        for i in range( word_length, len( ints ) ):
            index = 0
            factor = 1
            skip = False
            for j in range( word_length ):
                assert 0 < i-j < len( ints )
                letter = ints[i-j]
                if letter < 0:
                    skip = True
                    break
                index += letter * factor
                factor *= radix
            if skip: 
                continue        
            else:
                counts[ index ] += 1
                total += 1
        # Write ints separated by tabs
        print '\t'.join( [ str( total ) ] + map( str, counts ) )

if __name__ == "__main__": main()
