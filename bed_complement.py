#!/usr/bin/env python

"""
Complement the regions of a bed file

usage: %prog bed_file chrom_length_file
"""
import sys
from warnings import warn
from bx.bitset import *
import cookbook.doc_optparse

def read_len( f ):
    """Read a 'LEN' file and return a mapping from chromosome to length"""
    mapping = dict()
    for line in f:
        fields = line.split()
        mapping[ fields[0] ] = int( fields[1] )
    return mapping

options, args = cookbook.doc_optparse.parse( __doc__ )
try:
    in_fname, len_fname = args
except:
    cookbook.doc_optparse.exit()

bitsets = binned_bitsets_from_file( open( in_fname ) )

lens = read_len( open( len_fname ) )

# chrom_col, start_col, end_col = 0, 1, 2
# 
# last_chrom = None
# last_bitset = None
# bitsets = dict() 
# for line in open( in_fname ):
#     fields = line.split()
#     chrom = fields[chrom_col]
#     if chrom != last_chrom:
#         if chrom not in bitsets:
#             bitsets[chrom] = BitSet( lens[chrom] ) 
#         last_chrom = chrom
#         last_bitset = bitsets[chrom]
#     start, end = int( fields[start_col] ), int( fields[end_col] )
#     if start > end: warn( "Interval start after end!" )
#     last_bitset.set_range( start, end-start )

for chrom in lens:
    if chrom in bitsets:
        bits = bitsets[chrom]
        bits.invert()
        len = lens[chrom]
        end = 0
        while 1:
            start = bits.next_set( end )
            if start == bits.size: break
            end = bits.next_clear( start )
            if end > len: end = len
            print "%s\t%d\t%d" % ( chrom, start, end )
            if end == len: break
    else:
        print "%s\t%d\t%d" % ( chrom, 0, lens[chrom] )
