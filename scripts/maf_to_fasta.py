#!/usr/bin/env python.

"""
Read a maf and print the text as a fasta file.

usage: %prog < maf > fasta
"""

from __future__ import division

import textwrap
import sys
from bx.align import maf

def __main__():

    maf_reader = maf.Reader( sys.stdin )

    # Confusing since maf_to_concat_fasta takes names.

    # if len( sys.argv ) > 1:
    #     comps = map( int, sys.argv[1:] )
    # else:
    #     comps = None    

    comps = None

    for i, m in enumerate( maf_reader ):
        if comps: l = [ m.components[i] for i in comps ]
        else: l = m.components
        for c in l:
            print ">%s:%d-%d" % ( c.src, c.start, c.end )
            print c.text
            #print_n( c.text, 50 )

def print_n( s, n, f = sys.stdout ):
    p = 0
    while p < len( s ):
        print >> f, s[p:min(p+n,len(s))]
        p += n

if __name__ == "__main__": __main__()
