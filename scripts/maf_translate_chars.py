#!/usr/bin/env python2.4

from __future__ import division

import psyco_full

import sys

import sys
from bx.align import maf
from optparse import OptionParser
import string

table = string.maketrans( "#=X@", "-***")

def main():

    maf_reader = maf.Reader( sys.stdin )
    maf_writer = maf.Writer( sys.stdout )

    for m in maf_reader:
        for c in m.components:
            c.text = c.text.translate( table )
        maf_writer.write( m )    
    
    maf_writer.close()
    
if __name__ == "__main__": 
    main()
