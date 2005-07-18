#!/usr/bin/env python2.3

import psyco_full

import sys

import sys
from bx import align

def __main__():

    maf_reader = align.maf.Reader( sys.stdin )
    maf_writer = align.maf.Writer( sys.stdout )

    for m in maf_reader:
    
        align.shuffle_columns( m )

        maf_writer.write( m )

if __name__ == "__main__": __main__()
