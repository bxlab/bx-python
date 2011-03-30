"""
Concatenate sets of intervals. 

Preserves format of the first input -- it is possible to concat two files that
have different column orders. Of course, the meta-data of the second will be
lost (and filled with a "."). If all of the files (GenomicInteralReaders) are
the same format, sameformat=True will preserve all columns of the first input,
cuts extra columns on subsequent input, and pads missing columns. If
sameformat=False then extra columns are filled with ".".
"""

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def concat(readers, comments=True, header=True, sameformat=True):
    # Save columns from the first input
    chrom_col = readers[0].chrom_col
    start_col = readers[0].start_col
    end_col = readers[0].end_col
    strand_col = readers[0].strand_col
    nfields = None
    firstdataset = True
    output = False
    for intervals in readers:
        for interval in intervals:
            if isinstance(interval, GenomicInterval):
                if not nfields: nfields = interval.nfields
                out_interval = interval.copy()
                if sameformat or firstdataset:
                    # everything except the first input has to be
                    # trimmed or padded to match the first input
                    if len(out_interval.fields) > nfields:
                        out_interval.fields = out_interval.fields[0:nfields]
                        while len(out_interval.fields) < nfields:
                            out_interval.fields.append(".")
                    output = True
                    yield out_interval
                else:
                    chrom = out_interval.chrom
                    start = out_interval.start
                    end = out_interval.end
                    strand = out_interval.strand
                    out_interval.fields = ["." for col in range(nfields)]  
                    out_interval.fields[chrom_col] = chrom
                    out_interval.fields[start_col] = str(start)
                    out_interval.fields[end_col] = str(end)
                    # Strand is optional, might not exist in output
                    if strand_col < len( out_interval.fields ):
                        out_interval.fields[strand_col] = strand
                    yield out_interval
            elif isinstance(interval, Header) and header:
                yield interval
            elif isinstance(interval, Comment) and comments:
                yield interval
        if output and firstdataset: firstdataset = False
