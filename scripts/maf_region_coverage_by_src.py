#!/usr/bin/env python

"""
Reads a list of intervals and a set of indexed mafs. For each interval print
the amount covered by each species other than the reference.

usage: %prog maf_files  [options] < interval_file
   -s, --src=s:      Use this src for all intervals
   -p, --prefix=p:   Prepend this to each src before lookup
"""

import sys
from collections import defaultdict

import bx.align.maf
from bx.cookbook import doc_optparse


def __main__():

    # Parse Command Line

    options, args = doc_optparse.parse(__doc__)

    try:
        maf_files = args
        if options.prefix:
            prefix = options.prefix
        else:
            prefix = None
    except Exception:
        doc_optparse.exit()

    # Open indexed access to mafs
    indexes = [bx.align.maf.Indexed(maf_file, maf_file + ".index") for maf_file in maf_files]

    # Iterate over input ranges

    for line in sys.stdin:
        fields = line.split()
        src, start, end = fields[0], int(fields[1]), int(fields[2])
        if prefix:
            src = prefix + src

        total_length = end - start

        # Find overlap with reference component
        blocks = []
        for index in indexes:
            blocks += index.get(src, start, end)

        coverage = defaultdict(int)
        for block in blocks:
            overlap_start = max(start, block.components[0].start)
            overlap_end = min(end, block.components[0].end)
            length = overlap_end - overlap_start
            assert length > 0
            for c in block.components[1:]:
                species = c.src.split('.')[0]
                coverage[species] += length

        print(line, end=' ')
        for key, value in coverage.items():
            print("   ", key.ljust(10), "%0.2f" % (value / total_length))


if __name__ == "__main__":
    __main__()
