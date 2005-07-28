#!/usr/bin/env python

"""
Find regions of one bed file that overlap regions in a second bed file
usage: %prog bed_file_1 bed_file_2
"""
import sys
from warnings import warn
from bx.bitset import BinnedBitSet

in_fname, in2_fname = sys.argv[1:]

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
    ns = bitsets[fields[0]].next_set( start )
    if ns is not None and ns < end: 
        print line,
