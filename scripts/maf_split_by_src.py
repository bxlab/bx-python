#!/usr/bin/env python2.3

"""
Read a MAF from stdin and break into several mafs based on the source of 
each block. If the `component` option is provided then only that component
will be used to determine the new file for each block, otherwise the src
for *all* components will be used.

TODO: Should be able to specify component by species/prefix?

usage: %prog [options] < maf
    -o, --outprefix: prepend this to the name of each generate maf
    -c, --component: use only this component (by index!) to split
"""

import sys, string
import bx.align.maf
from optparse import OptionParser

import psyco_full

INF="inf"

def __main__():

    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "-o", "--outprefix", action="store", default="" )
    parser.add_option( "-c", "--component", action="store", default=None )
    ( options, args ) = parser.parse_args()

    out_prefix = options.outprefix
    comp = options.component
    if comp is not None:
       comp = int( comp )

    maf_reader = bx.align.maf.Reader( sys.stdin )

    writers = {}

    for m in maf_reader:
       
        if comp is None: 
            writer_key = string.join( [ c.src for c in m.components ], '_' )
        else:
            writer_key = m.components[ comp ].src

        if not writers.has_key( writer_key ):
            writer = bx.align.maf.Writer( file( "%s%s.maf" % ( out_prefix, writer_key ), "w" ) )
            writers[ writer_key ] = writer
        else:
            writer = writers[ writer_key ] 

        writer.write( m )

    for key in writers:
        writers[ key ].close()

if __name__ == "__main__": __main__()
