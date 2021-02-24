#!/usr/bin/env python

"""
Given two bed files print the number of bases covered 1) by both, 2) only by
the first, and 3) only by the second.

usage: %prog bed_file_1 bed_file_2
"""

from bx.bitset_builders import binned_bitsets_from_file
from bx.cookbook import doc_optparse


def coverage(bitsets):
    total = 0
    for chrom in bitsets:
        total += bitsets[chrom].count_range(0, bitsets[chrom].size)
    return total


options, args = doc_optparse.parse(__doc__)
try:
    in_fname, in2_fname = args
except ValueError:
    doc_optparse.exit()

bits1 = binned_bitsets_from_file(open(in_fname))
bits2 = binned_bitsets_from_file(open(in2_fname))

bits1_covered = coverage(bits1)
bits2_covered = coverage(bits2)

bitsets = dict()

for key in bits1:
    if key in bits2:
        bits1[key].iand(bits2[key])
        bitsets[key] = bits1[key]

both_covered = coverage(bitsets)

print("in both:  \t%d" % both_covered)
print("only in %s:\t%d" % (in_fname, bits1_covered - both_covered))
print("only in %s:\t%d" % (in2_fname, bits2_covered - both_covered))
