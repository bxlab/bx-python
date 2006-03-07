#!/usr/bin/env python2.3

"""
Read a maf and print the text as a fasta file, concatenating blocks

usage %prog species1,species2 < maf_file
"""

from __future__ import division

import textwrap
import sys
from bx.align import maf

def __main__():
    species = sys.argv[1].split(',')
    texts = {}
    for s in species: texts[s] = []
    maf_reader = maf.Reader( sys.stdin )
    for m in maf_reader:
        for s in species:
            c = m.get_component_by_src_start( s ) 
            if c: texts[s].append( c.text )
            else: texts[s].append( "-" * m.text_len )
    for s in species:
        print ">" + s
        print_n( "".join( texts[s] ), 50 )

def print_n( s, n, f = sys.stdout ):
    p = 0
    while p < len( s ):
        print >> f, s[p:min(p+n,len(s))]
        p += n

if __name__ == "__main__": __main__()
