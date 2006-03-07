#!/usr/bin/env python

"""
Find regions of first bed file that overlap regions in a second bed file

usage: %prog bed_file_1 bed_file_2
    -m, --mincols=N: Require this much overlap (default 1bp)
    -u, --upstream_pad=N: upstream interval padding (default 0bp)
    -d, --downstream_pad=N: downstream interval padding (default 0bp)
    -v, --reverse: Print regions that DO NOT overlap
"""

import sys
from warnings import warn

from bx.bitset import *
from bx.bitset_builders import *

import cookbook.doc_optparse

import pkg_resources
pkg_resources.require( "bx-python" )

mincols = 1
upstream_pad = 0
downstream_pad = 0

options, args = cookbook.doc_optparse.parse( __doc__ )
try:
    if options.mincols: mincols = int( options.mincols )
    if options.upstream_pad: upstream_pad = int( options.upstream_pad )
    if options.downstream_pad: downstream_pad = int( options.downstream_pad )
    reverse = bool( options.reverse )
    in_fname, in2_fname = args
except:
    cookbook.doc_optparse.exit()

# Read first bed into some bitsets

bitsets = binned_bitsets_from_file( open( in2_fname ) )

# Read second BED and intersect

for line in open( in_fname ):
    fields = line.split()
    if fields[0] not in bitsets: continue
    start, end = int( fields[1] ), int( fields[2] )
    if start > end: warn( "Bed interval start after end!" )
    if bitsets[fields[0]].count_range( start, end-start ) >= mincols:
        if not reverse: print line,
    elif reverse:
        print line,
