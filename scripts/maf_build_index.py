#!/usr/bin/env python2.3

"""
Build an index file for a set of MAF alignment blocks.

If index_file is not provided maf_file.index is used.

usage: %prog maf_file index_file
"""

import psyco_full

import cookbook.doc_optparse

from bx import interval_index_file
import sys

import bx.align.maf

def main():

    # Parse command line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        maf_file = sys.argv[1]
        if len( sys.argv ) > 2: index_file = sys.argv[2]
        else: index_file = maf_file + ".index" 
    except:
        cookbook.doc_optparse.exit()

    maf_reader = bx.align.maf.Reader( open( maf_file ) )

    indexes = interval_index_file.Indexes()

    # Need to be a bit tricky in our iteration here to get the 'tells' right
    while 1:
        pos = maf_reader.file.tell()
        block = maf_reader.next()
        if block is None: break
        for c in block.components:
            indexes.add( c.src, c.forward_strand_start, c.forward_strand_end, pos )

    out = open( index_file, 'w' )
    indexes.write( out )
    out.close()

if __name__ == "__main__": main()
