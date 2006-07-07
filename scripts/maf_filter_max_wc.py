#!/usr/bin/env python2.4

from __future__ import division

import psyco_full

import sys

import ranges, sys
from bx.align import maf
from optparse import OptionParser

def main():

    min_good = int( sys.argv[1] )
    min_species = int( sys.argv[2] )

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    for m in maf_reader:
        good = 0
        for col in m.column_iter():
            if col.count( '*' ) <= min_species:
                good += 1   
        if good >= min_good: 
            maf_writer.write( m )

if __name__ == "__main__": 
    main()
