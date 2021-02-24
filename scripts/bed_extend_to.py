#!/usr/bin/env python

"""
Read BED file and extend each record to the specified minimum length. If chromosome
size information is provided trim extended intervals.

usage: %prog amount [ chrom_file ] < bed_file
"""

import sys

from bx.intervals.io import GenomicIntervalReader

length = int(sys.argv[1])
chrom_len = None
if len(sys.argv) > 2:
    chrom_len = {fields[0]: int(fields[1]) for fields in map(str.split, open(sys.argv[2]))}

for interval in GenomicIntervalReader(sys.stdin):
    if interval.end - interval.start < length:
        start = interval.start
        end = interval.end
        # Extend in positive direction on strand
        if interval.strand == "+":
            end = start + length
        else:
            start = end - length
        # Trim
        if start < 0:
            start = 0
        if chrom_len and end > chrom_len[interval.chrom]:
            end = chrom_len[interval.chrom]
        # Set new start and end
        interval.start = start
        interval.end = end
    # Output possibly adjusted interval
    print(interval)
