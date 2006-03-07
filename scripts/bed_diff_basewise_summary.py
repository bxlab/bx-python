#!/usr/bin/env python

"""
Find regions of first bed file that overlap regions in a second bed file

usage: %prog bed_file_1 bed_file_2
"""
import sys
from warnings import warn
from bx.bitset import BinnedBitSet
from bx.bitset_builders import *
import cookbook.doc_optparse

def coverage( bitsets ):
    total = 0
    for chrom in bitsets:
        total += bitsets[chrom].count_range( 0, bitsets[chrom].size )
    return total    

options, args = cookbook.doc_optparse.parse( __doc__ )
try:
    in_fname, in2_fname = args
except:
    cookbook.doc_optparse.exit()

bits1 = binned_bitsets_from_file( open( in_fname ) )
bits2 = binned_bitsets_from_file( open( in2_fname ) )

bits1_covered = coverage( bits1 )
bits2_covered = coverage( bits2 )

bitsets = dict()

for key in bits1:
    if key in bits2:
        bits1[key].iand( bits2[key] )
        bitsets[key] = bits1[key]

both_covered = coverage( bitsets )

print "in both:  \t%d" % both_covered
print "only in %s:\t%d" % ( in_fname, bits1_covered - both_covered )
print "only in %s:\t%d" % ( in2_fname, bits2_covered - both_covered )
