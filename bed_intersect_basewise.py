#!/usr/bin/env python

"""
Find regions of first bed file that overlap regions in a second bed file

usage: %prog bed_file_1 bed_file_2
"""
import sys
from warnings import warn
from bx.bitset import BinnedBitSet
import cookbook.doc_optparse

def read_bed( f ):
    last_chrom = None
    last_bitset = None
    bitsets = dict() 
    for line in f:
        fields = line.split()
        if fields[0] != last_chrom:
            if fields[0] not in bitsets:
                bitsets[fields[0]] = BinnedBitSet() 
            last_chrom = fields[0]
            last_bitset = bitsets[fields[0]]
        start, end = int( fields[1] ), int( fields[2] )
        if start > end: warn( "Bed interval start after end!" )
        last_bitset.set_range( start, end-start )
    return bitsets   

options, args = cookbook.doc_optparse.parse( __doc__ )
try:
    in_fname, in2_fname = args
except:
    cookbook.doc_optparse.exit()

bits1 = read_bed( open( in_fname ) )
bits2 = read_bed( open( in2_fname ) )

bitsets = dict()

for key in bits1:
    if key in bits2:
        bits1[key].iand( bits2[key] )
        bitsets[key] = bits1[key]

for chrom in bitsets:
    bits = bitsets[chrom]
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print "%s\t%d\t%d" % ( chrom, start, end )
