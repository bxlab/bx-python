#!/usr/bin/env python

"""
For every column that occurs in a multiple alignment print the column
and the number of times it occurs (one column/count per line, tab
separated), sorted by count descending.

Note: all blocks must have exactly the same number of species.

usage: %prog < maf > column_counts
"""

import bx.align.maf
import sys

from itertools import *

counts = {}

nspecies = None

for block in bx.align.maf.Reader( sys.stdin ):
    # Ensure all blocks have the same number of rows
    if nspecies: assert len( block.components ) == nspecies
    else: nspecies = len( block.components )
    # Increment count for each column
    for col in izip( * [ iter( comp.text.upper() ) for comp in block.components ] ):
        try: counts[ col ] += 1
        except: counts[ col ] = 1

counts = [ ( value, key ) for key, value in counts.iteritems() ]
counts.sort()
counts.reverse()

# print len( counts )

for count, col in counts:
    print "".join(col), count
