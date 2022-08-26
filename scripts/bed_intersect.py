#!/usr/bin/env python

"""
Find regions of first bed file that overlap regions in a second bed file. The
output preserves all fields from the input.

NOTE: -u and -d options are currently not functional!

usage: %prog bed_file_1 bed_file_2
    -m, --mincols=N: Require this much overlap (default 1bp)
    -u, --upstream_pad=N: upstream interval padding (default 0bp)
    -d, --downstream_pad=N: downstream interval padding (default 0bp)
    -v, --reverse: Print regions that DO NOT overlap
    -b, --booleans: Just print '1' if interval overlaps or '0' otherwise
"""

from warnings import warn

from bx.bitset_builders import binned_bitsets_from_file
from bx.cookbook import doc_optparse

mincols = 1
upstream_pad = 0
downstream_pad = 0

options, args = doc_optparse.parse(__doc__)
try:
    if options.mincols:
        mincols = int(options.mincols)
    if options.upstream_pad:
        upstream_pad = int(options.upstream_pad)
    if options.downstream_pad:
        downstream_pad = int(options.downstream_pad)
    reverse = bool(options.reverse)
    booleans = bool(options.booleans)
    in_fname, in2_fname = args
except Exception:
    doc_optparse.exit()

# Read first bed into some bitsets

bitsets = binned_bitsets_from_file(open(in2_fname))

# Read second BED and intersect

for line in open(in_fname):
    if line.startswith("#") or line.isspace():
        continue
    fields = line.split()
    start, end = int(fields[1]), int(fields[2])
    if start > end:
        warn("Bed interval start after end!")
    if fields[0] in bitsets and bitsets[fields[0]].count_range(start, end - start) >= mincols:
        if booleans:
            if reverse:
                print(0)
            else:
                print(1)
        elif not reverse:
            print(line, end=" ")
    else:
        if booleans:
            if reverse:
                print(1)
            else:
                print(0)
        elif reverse:
            print(line, end=" ")
