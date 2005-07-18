#!/usr/bin/env python2.3

"""
Read a maf from stdin and print the chromosome number or each alignment

usage: %prog refindex [options]
"""

from __future__ import division

import sys
import cookbook.doc_optparse
from bx.align import maf
from optparse import OptionParser

def __main__():

    # Parse command line arguments
    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        refindex = int( args[0] )
    except:
        cookbook.doc_optparse.exit()

    maf_reader = maf.Reader( sys.stdin )

    for m in maf_reader: 
		c = m.components[ refindex ].src
		print c[ c.rfind( "chr" ) + 3 : ]

if __name__ == "__main__": __main__()
