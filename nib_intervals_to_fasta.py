#!/usr/bin/env python2.3

"""
Read a set of ranges and a nib file, print portions of nib overlapping
those ranges to stdout

usage: %prog range_file nib_file
"""

import cookbook.doc_optparse
from bx import seq.nib
import string
import sys

def __main__():

    options, args = cookbook.doc_optparse.parse( __doc__ )

    try:
        range_file = file( args[0] )
        nib_file = file( args[1] )
    except:
        cookbook.doc_optparse.exit()

    nib = seq.nib.NibFile( nib_file )

    for line in range_file: 
        fields = line.split()
        start, end = int( fields[0] ), int( fields[1] ) 
        print ">", start, end 
        print_wrapped( nib.get( start, end - start ) )

def print_wrapped( s ):
    l = len( s )        
    c = 0
    while c < l:
        b = min( c + 50, l )
        print s[c:b]
        c = b

if __name__ == "__main__": __main__()
