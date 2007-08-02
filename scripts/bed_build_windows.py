#!/usr/bin/env python2.4

"""
Build windows of length `window_size` over the sequences defined by 
`len_file` excluding regions in `gap_file`.

After removing the gaps, windows of exactly `window_size` units will be
placed in the remaining regions, with the extra space evenly placed
between the windows.

`len_file` is LEN format (name length) and `gap_file is BED (name start stop).

usage: %prog len_file gap_file window_size
"""

import sys
from bx.bitset_builders import binned_bitsets_from_file
import random

def main():
    region_fname, exclude_fname, window_size = sys.argv[1], sys.argv[2], int( sys.argv[3] ) 
    exclude_bitsets = binned_bitsets_from_file( open( exclude_fname ) )
    for line in open( region_fname ):
        fields = line.split()
        chr, start, end = fields[0], 0, int( fields[1] )
        if chr not in exclude_bitsets:
            do_windows( chr, start, end, window_size )
        else:
            bits = exclude_bitsets[chr]
            assert end < bits.size
            e = 0
            while 1:
                s = bits.next_clear( e )
                if s > end: break
                e = bits.next_set( s )
                do_windows( chr, s, min( e, end ), window_size )

def do_windows( chr, start, end, window_size ):
    length = end - start
    window_count = length // window_size
    if window_count == 0:
        return
    lost = length % window_size
    skip_amount = lost // window_count
    ## skip_amounts = [0] * ( window_count + 1 )
    ## for i in range( 0, lost ): skip_amounts[ random.randrange( 0, window_count + 1 ) ] += 1
    s = 0
    for i in range( 0, window_count ):
        s += skip_amount
        print chr, start + s, start + s + window_size        
        s += window_size
        
if __name__ == "__main__":
    main()