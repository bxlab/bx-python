#!/usr/bin/env python2.4

"""
Takes a list of maf filenames on the command line and prints a comma separated
list of the species that occur in all of the mafs. 

usage %prog maf1 maf2 ...
"""

import operator
import sys
import bx.align.maf

files = sys.argv[1:]
sets = []

for file in files:
    sys.stderr.write( "." )
    s = set()
    for block in bx.align.maf.Reader( open( file ) ):
        for comp in block.components:
            s.add( comp.src.split( '.' )[0] )
    sets.append( s )

inter = reduce( operator.and_, sets )
print ",".join( inter )
