#!/usr/bin/env python

"""
Build an index file for a set of MAF alignment blocks.

If index_file is not provided maf_file.index is used.

usage: %prog maf_file index_file
    -s, --species=a,b,c: only index the position of the block in the listed species
"""

import psyco_full

from bx.cookbook import doc_optparse

from bx import interval_index_file
import sys

import bx.align.maf

def main():

    # Parse command line

    options, args = doc_optparse.parse( __doc__ )

    try:
        maf_file = args[0]
        if len( args ) > 1: index_file = args[1]
        else: index_file = maf_file + ".index" 
        if options.species:
            species = options.species.split( "," )
        else:
            species = None
    except:
        doc_optparse.exit()

    maf_reader = bx.align.maf.Reader( open( maf_file ) )

    indexes = interval_index_file.Indexes()

    # Need to be a bit tricky in our iteration here to get the 'tells' right
    while 1:
        pos = maf_reader.file.tell()
        block = maf_reader.next()
        if block is None: break
        for c in block.components:
            if species is not None and c.src.split('.')[0] not in species:
                continue
            indexes.add( c.src, c.forward_strand_start, c.forward_strand_end, pos )

    out = open( index_file, 'w' )
    indexes.write( out )
    out.close()

if __name__ == "__main__": main()
