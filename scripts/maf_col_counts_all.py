#!/usr/bin/env python

"""
For every column that occurs in a multiple alignment print the column
and the number of times it occurs (one column/count per line, tab
separated), sorted by count descending.

This version allows special handling of the 'wildcard' symbol in alignments.

Note: all blocks must have exactly the same number of species.

usage: %prog [options] < maf > column_counts
    -w, --wildcard: include wildcards
    -m, --maxwildcards=N: only allow N missing species
"""

import bx.align.maf
import sys
from bx.cookbook import doc_optparse, cross_lists

from itertools import *

counts = {}

nspecies = None

for block in bx.align.maf.Reader( sys.stdin ):
    # Ensure all blocks have the same number of rows
    if nspecies: assert len( block.components ) == nspecies
    else: nspecies = len( block.components )
    # Increment count for each column
    for col in izip( * [ iter( comp.text.upper() ) for comp in block.components ] ):
        col = ''.join( col )
        try: counts[ col ] += 1
        except: counts[ col ] = 1

# counts = [ ( value, key ) for key, value in counts.iteritems() ]
# counts.sort()
# counts.reverse()

## for count, col in counts:
##     print "".join(col), count

options, args = doc_optparse.parse( __doc__ )

wildcard = False
if options.wildcard: 
    wildcard = True
    max_wildcard = nspecies - 1
if options.maxwildcards:
    wildcard = True
    max_wildcard = int( options.maxwildcards ) 

nucs = "ACGT-"
if wildcard:
    nucs += "*"

for col in cross_lists( *( [ nucs ] * nspecies ) ):
    col = ''.join( col )
    if wildcard and col.count( "*" ) > max_wildcard:
        continue
    if col.count( "-" ) == nspecies:
        continue
    print col, counts.get( col, 0 )
