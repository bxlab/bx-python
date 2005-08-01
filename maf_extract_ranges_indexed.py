#!/usr/bin/env python2.3

"""
Reads a list of intervals and a maf. Produces a new maf containing the
blocks or parts of blocks in the original that overlapped the intervals.

If index_file is not provided maf_file.index is used.

NOTE: If two intervals overlap the same block it will be written twice. With
      non-overlapping intervals and --chop this is never a problem. 

usage: %prog maf_file index_file [options] < interval_file
   -m, --mincols=10: Minimum length (columns) required for alignment to be output
   -c, --chop:       Should blocks be chopped to only portion overlapping (no by default)
   -s, --src=s:      Use this src for all intervals
   -p, --prefix=p:   Prepend this to each src before lookup
   -d, --dir=d:      Write each interval as a separate file in this directory
   -S, --strand:     Strand is included as an additional column, and the blocks are reverse complemented so that they are always on the plus strand w/r/t the src species.
"""

import psyco_full

import cookbook.doc_optparse

import bx.align.maf
from bx import misc
import os
import sys

def __main__():

    # Parse Command Line

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        maf_files = args
        if options.mincols: mincols = int( options.mincols )
        else: mincols = 10
        if options.src: fixed_src = options.src
        else: fixed_src = None
        if options.prefix: prefix = options.prefix
        else: prefix = None
        if options.dir: dir = options.dir
        else: dir = None
        chop = bool( options.chop )
        do_strand = bool( options.strand )
    except:
        cookbook.doc_optparse.exit()

    # Open indexed access to mafs
    index = bx.align.maf.MultiIndexed( maf_files )

    # Start MAF on stdout

    if dir is None: 
        out = bx.align.maf.Writer( sys.stdout )

    # Iterate over input ranges 

    for line in sys.stdin:
        strand = "+"
        fields = line.split()
        if fixed_src:
            src, start, end = fixed_src, int( fields[0] ), int( fields[1] )
            if do_strand: strand = fields[2]
        else:
            src, start, end = fields[0], int( fields[1] ), int( fields[2] )
            if do_strand: strand = fields[3]
        if prefix: src = prefix + src
        # Find overlap with reference component
        blocks = index.get( src, start, end )
        # Open file if needed
        if dir:
            out = bx.align.maf.Writer( open( os.path.join( dir, "%s:%09d-%09d.maf" % ( src, start, end ) ), 'w' ) )
        # Write each intersecting block
        if chop:
            for block in blocks: 
                ref = block.get_component_by_src( src )
                # If the reference component is on the '-' strand we should complement the interval
                if ref.strand == '-':
                    slice_start = max( ref.src_size - end, ref.start )
                    slice_end = max( ref.src_size - start, ref.end )
                else:
                    slice_start = max( start, ref.start )
                    slice_end = min( end, ref.end )
                sliced = block.slice_by_component( ref, slice_start, slice_end ) 
                good = True
                for c in sliced.components: 
                    if c.size < 1: 
                        good = False
                if good and sliced.text_size > mincols:
                    if strand != ref.strand: sliced = sliced.reverse_complement()
                    out.write( sliced )
        else:
            for block in blocks:
                out.write( block )
        if dir:
            out.close()
         
    # Close output MAF

    out.close()

if __name__ == "__main__": __main__()
