#!/usr/bin/env python

"""
Pass through blocks from a maf file until a certain number of columns
have been passed.

usage: %prog -c cols < maf > maf
"""

import sys

from bx.align import maf
from optparse import OptionParser

def __main__():

    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "-c", "--cols",  action="store" )

    ( options, args ) = parser.parse_args()

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    if not options.cols: raise "Cols argument is required"
    cols = int( options.cols )

    count = 0

    for m in maf_reader:

        maf_writer.write( m )

        count += m.text_size

        if count >= cols: return        

if __name__ == "__main__": __main__()
