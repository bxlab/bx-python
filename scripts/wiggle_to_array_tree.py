#!/usr/bin/env python

"""
Read data in UCSC wiggle format and write it to an "array tree" file.

usage: %prog array_length output.array_tree < input.wig
"""

from __future__ import division

import sys

from bx.arrays.array_tree import ArrayTree
from bx.arrays.wiggle import IntervalReader

def main():
   
    max = int( sys.argv[1] )
    out_fname = sys.argv[2]
    
    scores = ArrayTree( max, 1000 )

    # Fill array from wiggle
    for i, ( chrom, start, end, _, val ) in enumerate( IntervalReader( sys.stdin ) ):
        scores.set_range( start, end, val )
        if i % 100000 == 0: 
            print i, "scores processed -- last value:", scores[end-1]

    scores.root.build_summary()
    
    f = open( out_fname, "w" )
    scores.to_file( f )
    f.close()

if __name__ == "__main__": 
    main()
