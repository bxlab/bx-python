#!/usr/bin/env python2.3

"""
For each block in a maf file (read from stdin) write a sequence of ints 
corresponding to the columns of the block after applying the provided mapping.

The 'correct' number of species is determined by the mapping file, blocks not having
this number of species will be ignored.

usage: %prog mapping_file
"""

from __future__ import division

import psyco_full

from numarray import *

import bx.align.maf
from bx import alphabet
from bx import seq_numarray
import string
import sys

def main():

    if len( sys.argv ) > 1:
        alpha_map = alphabet.Mapping( file( sys.argv[1] ) )
    else:
        alpha_map = None

    for maf in bx.align.maf.Reader( sys.stdin ):
        # Translate alignment to ints
        int_seq = seq_numarray.DNA.translate_alignment( [ c.text for c in maf.components ] )
        # Apply mapping 
        if alpha_map:
            int_seq = alpha_map.translate( int_seq )
        # Write ints separated by spaces
        for i in int_seq: print i,
        print

if __name__ == "__main__": main()
