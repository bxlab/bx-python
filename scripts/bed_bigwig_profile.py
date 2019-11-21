#!/usr/bin/env python

"""
Create a site profile vector showing the average signal accumulated from a
bigwig file around the center of each interval from a BED file.

Output is the average signal value at that relative position across the
intervals.

usage: %prog bigwig_file.bw padding < bed_file.bed
"""

import sys

from numpy import (
    float64,
    floor,
    int32,
    isnan,
    savetxt,
    zeros
)

from bx.bbi.bigwig_file import BigWigFile
from bx.intervals.io import GenomicIntervalReader

bw = BigWigFile(open(sys.argv[1]))
padding = int(sys.argv[2])
totals = zeros(padding*2, dtype=float64)
valid = zeros(padding*2, dtype=int32)

for interval in GenomicIntervalReader(sys.stdin):
    center = floor((interval.start + interval.end) / 2)
    values = bw.get_as_array(interval.chrom, center - padding, center + padding)
    # Determine which positions had data and mask the rest for totalling
    invalid = isnan(values)
    values[invalid] = 0
    totals += values
    valid += (~ invalid)

savetxt(sys.stdout, totals/valid)
