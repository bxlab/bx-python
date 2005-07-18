#!/usr/bin/env python2.3

"""
Read a MAF from standard input and print counts of alignments, bases, or columns. 

usage: %prog [options]
   -l, --lnorm: divide (normalize) score by alignment text length
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
        lnorm = bool( options.lnorm )
    except:
        cookbook.doc_optparse.exit()

    maf_reader = maf.Reader( sys.stdin )

    for m in maf_reader: 
        if lnorm:
            print float( m.score ) / m.text_size
        else:
            print m.score

if __name__ == "__main__": __main__()
