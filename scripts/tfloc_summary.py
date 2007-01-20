#!/usr/bin/env python2.3

"""
Read TFLOC output from stdin and write out a summary in which the nth line 
contains the number of sites found in the nth alignment of the input.

TODO: This is very special case, should it be here?
"""

import sys

counts = dict()

max_index = -1

for line in sys.stdin:
    if line[0].isdigit():
        current_index = int( line )
        max_index = max( current_index, max_index )        
    elif line[0] == "'":
        try: counts[ current_index ] += 1
        except: counts[ current_index ] = 1
    else:
        raise "Invalid input line " + line

for i in range( max_index + 1 ):
    print counts.get( i, 0 )
