#!/usr/bin/env python

"""
Translate a maf file containing gap ambiguity characters as produced by 
'maf_tile_2.py' to a new file in which "#" (contiguous) is replaced by "-" and
all other types are replaces by "*".

TODO: This could be much more general, should just take the translation table
      from the command line.
      
usage: %prog < maf > maf
"""

from __future__ import division

import psyco_full

import sys

import sys
from bx.align import maf
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
