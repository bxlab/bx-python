#!/usr/bin/env python

"""
Filter each block in a maf file. Can filter blocks for a minimum number of
components (rows), a minimum length in columns, or an arbitrary python 
expression (which will be evaluated for each block with the variable 'm'
containing that block).

usage: %prog [options] < maf > maf
    --component_count=N: Minimum number of components (rows)
    --min_cols=N: Minimum number of columns
    -e, --expr=EXPR: Python expression that must evaulate to true
"""

import psyco_full

import sys

import sys
from bx.align import maf
from optparse import OptionParser

def __main__():

    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "--component_count", action="store", default=None, type="int", help="" )
    parser.add_option( "--min_cols", action="store", default=None, type="int", help="" )
    parser.add_option( "-e", "--expr", action="store", default=None )

    ( options, args ) = parser.parse_args()

    component_count = options.component_count
    min_cols = options.min_cols
    expr = options.expr

    # Compile expression for SPEED
    if expr: expr = compile( expr, '<expr arg>', 'eval' )

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    for m in maf_reader:

        if component_count and len( m.components ) != component_count: continue
        if min_cols and m.text_size < min_cols: continue
        if expr and not bool( eval( expr, { "m": m, "maf": m } ) ): continue

        maf_writer.write( m )

if __name__ == "__main__": __main__()
