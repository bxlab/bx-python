#!/usr/bin/env python

"""
Find continuous regions that are covered by the first bed file (`bed_file_1`)
but not by the second bed file (`bed_file_2`)

usage: %prog bed_file_1 bed_file_2
"""

import sys
from warnings import warn
from bx.bitset_builders import binned_bitsets_from_file
from bx.cookbook import doc_optparse

def print_bits_as_bed( bits ):
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print "%s\t%d\t%d" % ( chrom, start, end )

options, args = doc_optparse.parse( __doc__ )
try:
    in_fname, in2_fname = args
except:
    doc_optparse.exit()

# Read first bed into some bitsets

bitsets1 = binned_bitsets_from_file( open( in_fname ) )
bitsets2 = binned_bitsets_from_file( open( in2_fname ) )

for chrom in bitsets1:  
    if chrom not in bitsets1:
        continue
    bits1 = bitsets1[chrom]
    if chrom in bitsets2:
        bits2 = bitsets2[chrom]
        bits2.invert()
        bits1.iand( bits2 )
    print_bits_as_bed( bits1 )
    
