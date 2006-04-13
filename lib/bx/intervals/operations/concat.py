#!/usr/bin/env python
"""
Take an input of readers, and concatenate them.  Preserve format of
the first input.  In other words, it is possible to concat two files
that have different column orders.  Of course, the meta-data of the
second will be lost (and filled with a .).  If all of the files
(GenomicInteralReaders) are the same format, sameformat=True will
preserve all columns, but fails to pad extra columns.  If
sameformat=False then extra columns are filled with (.).
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def concat(readers, comments=True, header=True, sameformat=True):
    chrom_col = readers[0].chrom_col
    start_col = readers[0].start_col
    end_col = readers[0].end_col
    strand_col = readers[0].strand_col
    firsttime = True
    for intervals in readers:
        for interval in intervals:
            if type( interval ) is GenomicInterval:
                out_interval = interval.copy()
                if sameformat or firsttime:
                    yield out_interval
                    nfields = interval.nfields
                else:
                    chrom = out_interval.chrom
                    start = out_interval.start
                    end = out_interval.end
                    strand = out_interval.strand
                    out_interval.fields = ["."] * nfields
                    out_interval.fields[chrom_col] = chrom
                    out_interval.fields[start_col] = str(start)
                    out_interval.fields[end_col] = str(end)
                    out_interval.fields[strand_col] = strand
                    yield out_interval
            elif type( interval ) is Header and firsttime and header:
                yield interval
            elif type( interval ) is Comment and comments:
                yield interval
        firsttime=False
