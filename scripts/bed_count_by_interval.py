#!/usr/bin/env python
"""
For each interval in `bed1` count the number of intersecting regions in `bed2`.

usage: %prog bed1 bed2
"""

import sys

from bx.intervals import (
    Intersecter,
    Interval
)

bed1, bed2 = sys.argv[1:3]

ranges = {}
for line in open(bed2):
    fields = line.strip().split()
    chrom, start, end, = fields[0], int(fields[1]), int(fields[2])
    if chrom not in ranges:
        ranges[chrom] = Intersecter()
    ranges[chrom].add_interval(Interval(start, end))

for line in open(bed1):
    fields = line.strip().split()
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    other = " ".join(fields[3:])
    out = " ".join(fields[:3] + [other])
    if chrom in ranges:
        print(out, len(ranges[chrom].find(start, end)))
    else:
        print(out, 0)
