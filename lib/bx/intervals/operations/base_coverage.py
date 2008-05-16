"""
Determine the number of bases covered by a set of intervals.
"""

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *

def base_coverage(reader):
    bitsets = reader.binned_bitsets()
    coverage = 0
    for chrom in bitsets:
        try:
            coverage += bitsets[chrom].count_range(0, MAX_END)
        except IndexError, e:
            try:
                # This will work only if reader is a NiceReaderWrapper
                reader.skipped += 1
                # no reason to stuff an entire bad file into memmory
                if reader.skipped < 10:
                    reader.skipped_lines.append( ( reader.linenum, reader.current_line, str( e ) ) )
            except:
                pass
            continue
    return coverage
