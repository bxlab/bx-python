#!/usr/bin/env python

"""
Find regions of first bed that are not in second bed file (subtract
second from first)

usage: %prog bed_file_1 bed_file_2

"""
import sys
from warnings import warn
from bx.bitset_builders import binned_bitsets_from_file
import cookbook.doc_optparse

def read_len( f ):
    """Read a 'LEN' file and return a mapping from chromosome to length"""
    mapping = dict()
    for line in f:
        fields = line.split()
        mapping[ fields[0] ] = int( fields[1] )
    return mapping

def print_bits_as_bed( bits ):
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print "%s\t%d\t%d" % ( chrom, start, end )

options, args = cookbook.doc_optparse.parse( __doc__ )
try:
    in_fname, in2_fname, len_fname = args
except:
    cookbook.doc_optparse.exit()

lens = read_len( open( len_fname ) )

# Read first bed into some bitsets

bitsets1 = binned_bitsets_from_file( open( in_fname ), lens=lens )
bitsets2 = binned_bitsets_from_file( open( in2_fname ), lens=lens )

for chrom in bitsets1:  
    if chrom not in bitsets1:
        continue
    bits1 = bitsets1[chrom]
    if chrom in bitsets2:
        bits2 = bitsets2[chrom]
        bits2.invert()
        bits1.iand( bits2 )
    print_bits_as_bed( bits1 )
    
