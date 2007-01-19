#!/usr/bin/env python

"""
Find regions of first bed file that overlap regions in a second bed file. This
program performs a base-by-base intersection, so only runs of bases that are
covered in both of the inputs will be output.

usage: %prog bed_file_1 bed_file_2
"""

import sys
from warnings import warn
from bx.bitset import *
from bx.bitset_builders import *
from bx.cookbook import doc_optparse

options, args = doc_optparse.parse( __doc__ )
try:
    in_fname, in2_fname = args
except:
    doc_optparse.exit()

bits1 = binned_bitsets_from_file( open( in_fname ) )
bits2 = binned_bitsets_from_file( open( in2_fname ) )

bitsets = dict()

for key in bits1:
    if key in bits2:
        bits1[key].iand( bits2[key] )
        bitsets[key] = bits1[key]

for chrom in bitsets:
    bits = bitsets[chrom]
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print "%s\t%d\t%d" % ( chrom, start, end )
