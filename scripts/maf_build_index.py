#!/usr/bin/env python

"""
Build an index file for a set of MAF alignment blocks.

If index_file is not provided maf_file.index is used.

usage: %prog maf_file index_file
    -s, --species=a,b,c: only index the position of the block in the listed species
"""

import psyco_full

from bx.cookbook import doc_optparse

import sys
import os.path

from bx import interval_index_file
import bx.align.maf
from bx.misc.seekbzip2 import SeekableBzip2File

def main():

    # Parse command line

    options, args = doc_optparse.parse( __doc__ )

    try:
        maf_file = args[0]
        # If it appears to be a bz2 file, attempt to open with table
        if maf_file.endswith( ".bz2" ):
            table_file = maf_file + "t"
            if not os.path.exists( table_file ):
                doc_optparse.exit( "To index bz2 compressed files first "
                                   "create a bz2t file with bzip-table." )
            # Open with SeekableBzip2File so we have tell support
            maf_in = SeekableBzip2File( maf_file, table_file )
            # Strip .bz2 from the filename before adding ".index"
            maf_file = maf_file[:-4]
        elif maf_file.endswith( ".lzo" ):
            from bx.misc.seeklzop import SeekableLzopFile
            table_file = maf_file + "t"
            if not os.path.exists( table_file ):
                doc_optparse.exit( "To index lzo compressed files first "
                                   "create a lzot file with lzop_build_offset_table." )
            # Open with SeekableBzip2File so we have tell support
            maf_in = SeekableLzopFile( maf_file, table_file )
            # Strip .lzo from the filename before adding ".index"
            maf_file = maf_file[:-4]
        else:
            maf_in = open( maf_file )
        # Determine the name of the index file
        if len( args ) > 1: 
            index_file = args[1]
        else: 
            index_file = maf_file + ".index" 
        if options.species:
            species = options.species.split( "," )
        else:
            species = None
    except:
        doc_optparse.exception()

    maf_reader = bx.align.maf.Reader( maf_in )

    indexes = interval_index_file.Indexes()

    # Need to be a bit tricky in our iteration here to get the 'tells' right
    while 1:
        pos = maf_reader.file.tell()
        block = maf_reader.next()
        if block is None: break
        for c in block.components:
            if species is not None and c.src.split('.')[0] not in species:
                continue
            indexes.add( c.src, c.forward_strand_start, c.forward_strand_end, pos, max=c.src_size )

    out = open( index_file, 'w' )
    indexes.write( out )
    out.close()

if __name__ == "__main__": main()
