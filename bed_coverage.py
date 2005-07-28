#!/usr/bin/env python

"""
Print number of bases covered by intervals in a bed file

usage: %prog bed files ...
"""

import psyco_full
import sys
from bx.bitset import BinnedBitSet
from itertools import *

bed_filenames = sys.argv[1:]
if bed_filenames:
    input = chain( * imap( open, bed_filenames ) )
else:
    input = sys.stdin

last_chrom = None
last_bitset = None
bitsets = dict() 

for line in input:
        if line.startswith("#") or line.startswith("track"): continue
        fields = line.split()
        if fields[0] != last_chrom:
            if fields[0] not in bitsets:
                bitsets[fields[0]] = BinnedBitSet()
            last_chrom = fields[0]
            last_bitset = bitsets[fields[0]]
        start, end = int( fields[1] ), int( fields[2] )
        if start > end: print >>sys.stderr, "Bed interval start after end: " + line.strip()
        last_bitset.set_range( start, end-start )

total = 0
for chrom in bitsets:
    total += bitsets[chrom].count_range( 0, bitsets[chrom].size )

print total    
