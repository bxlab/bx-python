#!/usr/bin/env python

"""
Remove any blocks from a maf that overlap any of a set of intervals.

usage: %prog interval files... < maf
"""

import sys

import bx.align.maf
from bx import intervals
from bx.cookbook import doc_optparse


def __main__():
    options, args = doc_optparse.parse(__doc__)

    try:
        assert len(args) > 0
    except AssertionError:
        doc_optparse.exit()

    # Load Intervals

    intersector = intervals.Intersecter()

    for f in args:
        for line in open(f):
            if line.startswith("#") or line.isspace():
                continue
            fields = line.split()
            intersector.add_interval(intervals.Interval(int(fields[0]), int(fields[1])))

    # Start MAF on stdout

    out = bx.align.maf.Writer(sys.stdout)

    # Iterate over input MAF

    for maf in bx.align.maf.Reader(sys.stdin):
        # Find overlap with reference component
        intersections = intersector.find(maf.components[0].start, maf.components[0].end)
        # Write only if no overlap
        if len(intersections) == 0:
            out.write(maf)

    # Close output MAF

    out.close()


if __name__ == "__main__":
    __main__()
