#!/usr/bin/env python
"""
Take an input of a GenomicIntervalReader, merge overlapping regions,
and yield GenomicIntervals.
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

MAX_END = 512*1024*1024

# sorting could make this a less memory intensive operation(?)
def merge(interval, mincols=1):
    bitsets = interval.binned_bitsets()
    if interval.header:
        yield interval.header
    for chrom in bitsets:
        bitset = bitsets[chrom]
        output = ["."] * (max(interval.chrom_col, interval.start_col, interval.end_col) + 1)
        output[interval.chrom_col] = chrom
        for start, end in bits_set_in_range(bitset,0, MAX_END):
            output[interval.start_col] = str(start)
            output[interval.end_col] = str(end)
            yield output
