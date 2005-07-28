#!/usr/bin/env python

"""
Find regions of first bed file that overlap regions in a second bed file

usage: %prog bed_file_1 bed_file_2
    -m, --mincols=N: Require this much overlap 
"""
import sys
from warnings import warn
from bx.bitset import BinnedBitSet
import cookbook.doc_optparse

options, args = cookbook.doc_optparse.parse( __doc__ )
try:
    if options.mincols: mincols = int( options.mincols )
    else: mincols = 1
    in_fname, in2_fname = args
except:
    cookbook.doc_optparse.exit()


last_chrom = None
last_bitset = None
bitsets = dict() 

for line in open( in2_fname ):
    fields = line.split()
    if fields[0] != last_chrom:
        if fields[0] not in bitsets:
            bitsets[fields[0]] = BinnedBitSet() 
        last_chrom = fields[0]
        last_bitset = bitsets[fields[0]]
    start, end = int( fields[1] ), int( fields[2] )
    if start > end: warn( "Bed interval start after end!" )
    last_bitset.set_range( start, end-start )

# Read second BED and intersect

for line in open( in_fname ):
    fields = line.split()
    if fields[0] not in bitsets: continue
    start, end = int( fields[1] ), int( fields[2] )
    if start > end: warn( "Bed interval start after end!" )
    if bitsets[fields[0]].count_range( start, end-start ) >= mincols:
        print line,
