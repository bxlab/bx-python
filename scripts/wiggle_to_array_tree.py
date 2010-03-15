#!/usr/bin/env python

"""
Read data in UCSC wiggle format and write it to an "array tree" file.

usage: %prog array_length output.array_tree < input.wig
"""

from __future__ import division

import sys

from bx.arrays.array_tree import *
from bx.arrays.wiggle import WiggleReader

def main():
   
    sizes_fname = sys.argv[1]
    out_fname = sys.argv[2]
    
    sizes = {}
    for line in open( sizes_fname ):
        fields = line.split()
        sizes[ fields[0] ] = int( fields[1] )
    
    # Fill array from wiggle
    d = array_tree_dict_from_reader( WiggleReader( sys.stdin ), sizes )
    
    for value in d.itervalues():
        value.root.build_summary()
    
    f = open( out_fname, "w" )
    FileArrayTreeDict.dict_to_file( d, f )
    f.close()

if __name__ == "__main__": 
    main()
