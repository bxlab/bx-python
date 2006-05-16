#!/usr/bin/env python
"""
Complement a set of Intervals from reader.
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def complement(reader, lens):
    # load into bitsets
    bitsets = reader.binned_bitsets(upstream_pad = 0, downstream_pad = 0, lens = lens)

    # NOT them all
    for key, value in bitsets.items():
        value.invert()

    # Read remaining intervals and subtract
    for chrom in bitsets:
        bitset = bitsets[chrom]
        out_intervals = bits_set_in_range( bitset, 0, lens.get(chrom, 512*1024*1024))
        # Write the intervals
        for start, end in out_intervals:
            fields = ["."  for x in range(max(reader.chrom_col, reader.start_col, reader.end_col)+1)]
            fields[reader.chrom_col] = chrom
            fields[reader.start_col] = start
            fields[reader.end_col] = end
            new_interval = GenomicInterval(reader, fields, reader.chrom_col, reader.start_col, reader.end_col, reader.strand_col, "+")
            yield new_interval
