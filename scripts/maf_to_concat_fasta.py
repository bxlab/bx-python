#!/usr/bin/env python

"""
Read a maf and print the text as a fasta file, concatenating blocks. A 
specific subset of species can be chosen. 

usage %prog [options] species1,species2,... < maf_file > fasta_file
    --fill="expression": Insert this between blocks
    --wrap=columns: Wrap FASTA to this many columns
"""

from optparse import OptionParser

import textwrap
import sys
from bx.align import maf

def __main__():
    # Parse command line arguments

    parser = OptionParser()
    parser.add_option( "--fill", action="store", default=None, type="string", help="" )
    parser.add_option( "--wrap", action="store", default=None, type="int", help="" )
    parser.add_option( "--nowrap", action="store_true", default=False, dest="nowrap", help="" )

    ( options, args ) = parser.parse_args()

    species = []
    for arg in args: species.extend(arg.split(','))

    fill = ""
    if options.fill: fill = eval( options.fill )

    wrap = 50
    if   (options.wrap != None): wrap = options.wrap
    elif (options.nowrap):       wrap = 0

    # create the concatenated sequences

    texts = {}
    for s in species: texts[s] = []
    maf_reader = maf.Reader( sys.stdin )
    for m in maf_reader:
        for s in species:
            c = m.get_component_by_src_start( s ) 
            if c: texts[s].append( c.text )
            else: texts[s].append( "-" * m.text_size )
    for s in species:
        print ">" + s
        print_n( fill.join( texts[s] ), wrap )

def print_n( s, n, f = sys.stdout ):
    if (n <= 0):
        print >> f, s
    else:
        p = 0
        while p < len( s ):
            print >> f, s[p:min(p+n,len(s))]
            p += n

if __name__ == "__main__": __main__()
