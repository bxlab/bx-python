#!/usr/bin/env python2.4

"""
Read a MAF from standard input and print counts of alignments, bases, or columns. 

usage: %prog [options]
   -r, --recalculate: don't use the score from the maf, recalculate (using hox70 matrix)
   -l, --lnorm: divide (normalize) score by alignment text length
"""

from __future__ import division

import sys
import cookbook.doc_optparse
from bx.align import maf
from bx.align import score
from optparse import OptionParser

def main():

    # Parse command line arguments
    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        lnorm = bool( options.lnorm )
        recalculate = bool( options.recalculate )
    except:
        cookbook.doc_optparse.exit()

    hox70 = score.build_scoring_scheme( """  A    C    G    T
                                      91 -114  -31 -123
                                    -114  100 -125  -31
                                     -31 -125  100 -114
                                    -123  -31 -114   91 """, 400, 30, default=0 )

    maf_reader = maf.Reader( sys.stdin )

    for m in maf_reader: 
        if m.text_size == 0:
            print "NA"
            continue
        s = m.score
        # Recalculate?
        if recalculate:
            s = hox70.score_alignment( m )
        # Normalize?
        if lnorm:
            s = s / m.text_size
        # Print
        print s

if __name__ == "__main__": 
    main()
