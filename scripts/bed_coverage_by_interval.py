#!/usr/bin/env python

"""
For each interval in `bed1` print the fraction of bases covered by `bed2`.

usage: %prog bed1 bed2 [mask]
"""

import sys

from bx.bitset import BinnedBitSet
from bx.bitset_builders import binned_bitsets_from_file

bed1_fname, bed2_fname = sys.argv[1:3]

bitsets = binned_bitsets_from_file(open(bed2_fname))


def clone(bits):
    b = BinnedBitSet(bits.size)
    b.ior(bits)
    return b


if len(sys.argv) > 3:
    mask_fname = sys.argv[3]
    mask = binned_bitsets_from_file(open(mask_fname))
    new_bitsets = dict()
    for key in bitsets:
        if key in mask:
            b = clone(mask[key])
            b.invert()
            b.iand(bitsets[key])
            new_bitsets[key] = b
    bitsets = new_bitsets
else:
    mask = None

for line in open(bed1_fname):
    fields = line.split()
    chr, start, end = fields[0], int(fields[1]), int(fields[2])
    bases_covered = 0
    if chr in bitsets:
        bases_covered = bitsets[chr].count_range(start, end-start)
    length = end - start
    if mask and chr in mask:
        bases_masked = mask[chr].count_range(start, end-start)
        length -= bases_masked
    assert bases_covered <= length, f"{bases_covered!r}, {bases_masked!r}, {length!r}"
    if length == 0:
        print(0.0)
    else:
        print(bases_covered / length)
