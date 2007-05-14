#!/usr/bin/env python

"""
Read a set of ranges and a nib file, print portions of nib overlapping
those ranges to stdout

TODO: General sequence handling would be nice, as well as merging with
      'nib_intervals_to_fasta.py'.

usage: %prog nib_dir < range_file
"""

from bx.cookbook import doc_optparse
import bx.seq.nib
import string
import sys

def __main__():

    options, args = doc_optparse.parse( __doc__ )

    try:
        nib_dir = args[0] 
    except:
        doc_optparse.exit()

    nibs = {}

    for line in sys.stdin: 
        fields = line.split()
        chrom, start, end = fields[0], int( fields[1] ), int( fields[2] ) 
        print ">", chrom, start, end 
        if chrom in nibs:
            nib = nibs[chrom]
        else:
            nibs[chrom] = nib = bx.seq.nib.NibFile( file( "%s/%s.nib" % ( nib_dir, chrom ) ) )
        print_wrapped( nib.get( start, end - start ) )

def print_wrapped( s ):
    l = len( s )        
    c = 0
    while c < l:
        b = min( c + 50, l )
        print s[c:b]
        c = b

if __name__ == "__main__": __main__()
