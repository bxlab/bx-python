#!/usr/bin/env python

"""
Merge any overlapping regions of bed files. Bed files can be provided on the
command line or on stdin. Merged regions are always reported on the '+' 
strand, and any fields beyond chrom/start/stop are lost. 

usage: %prog bed files ...
"""

import psyco_full
import sys

from bx.bitset import *
from bx.bitset_builders import *
from itertools import *

bed_filenames = sys.argv[1:]
if bed_filenames:
    input = chain( * imap( open, bed_filenames ) )
else:
    input = sys.stdin

bitsets = binned_bitsets_from_bed_file( input )

for chrom in bitsets:
    bits = bitsets[chrom]
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print "%s\t%d\t%d" % ( chrom, start, end )
