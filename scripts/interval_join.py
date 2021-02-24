#!/usr/bin/env python

"""
Match up intersecting intervals from two files. This performs a "full join",
any pair of intervals with any basewise overlap will be printed side-by-side.

usage: %prog bed1 bed2
"""

import sys

import bx.intervals.intersection
import bx.intervals.io


def main():
    intersecters = {}

    # Read second set into intersecter
    for interval in bx.intervals.io.GenomicIntervalReader(open(sys.argv[2])):
        if interval.chrom not in intersecters:
            intersecters[interval.chrom] = bx.intervals.Intersecter()
        intersecters[interval.chrom].add_interval(interval)

    # Join with first set
    for interval in bx.intervals.io.GenomicIntervalReader(open(sys.argv[1])):
        if interval.chrom in intersecters:
            intersection = intersecters[interval.chrom].find(interval.start, interval.end)
            for interval2 in intersection:
                print("\t".join([str(interval), str(interval2)]))


if __name__ == "__main__":
    main()
