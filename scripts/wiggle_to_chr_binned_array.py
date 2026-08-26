#!/usr/bin/env python

"""
Writes compressed data from a wiggle file by chromosome.

usage: %prog score_file < wiggle_data
"""

import sys

import bx.wiggle
from bx.binned_array import BinnedArray


def main() -> None:
    scores: dict[str, BinnedArray] = {}
    with open(sys.argv[1]) as f:
        for i, (chrom, pos, val) in enumerate(bx.wiggle.Reader(f)):
            if chrom not in scores:
                scores[chrom] = BinnedArray()
            scores[chrom][pos] = val

            # Status
            if i % 10000 == 0:
                print(i, "scores processed")

    for chr, binned_array in scores.items():
        with open(chr, "w") as out:
            binned_array.to_file(out)


if __name__ == "__main__":
    main()
