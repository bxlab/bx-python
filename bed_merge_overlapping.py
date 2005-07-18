#!/usr/bin/env python

"""
Merge any overlapping regions of a bed file.
Uses bitset projection similar to 'featureBits'
"""

import sys
from bx.bitset import BinnedBitSet

in_fname, out_fname = sys.argv[1:]

last_chrom = None
last_bitset = None
bitsets = dict() 

for line in open( in_fname ):
    fields = line.split()
    if fields[0] != last_chrom:
        if fields[0] not in bitsets:
            bitsets[fields[0]] = BinnedBitSet()
        last_chrom = fields[0]
        last_bitset = bitsets[fields[0]]
    start, end = int( fields[1] ), int( fields[2] )
    assert start <= end, "Bed interval start must be less than end"
    last_bitset.set_range( start, end-start )

for chrom in bitsets:
    bits = bitsets[chrom]
    end = 0
    while 1:
        start = bits.next_set( end )
        if start is None: break
        end = bits.next_clear( start )
        print "%s\t%d\t%d" % ( chrom, start, end )
