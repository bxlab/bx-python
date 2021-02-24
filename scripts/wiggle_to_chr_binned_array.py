#!/usr/bin/env python

"""
Writes compressed data from a wiggle file by chromosome.

usage: %prog score_file < wiggle_data
"""

import sys

import bx.wiggle
from bx.binned_array import BinnedArray
from bx.cookbook import doc_optparse


def main():
    options, args = doc_optparse.parse(__doc__)

    scores = {}
    for i, (chrom, pos, val) in enumerate(bx.wiggle.Reader(open(sys.argv[1]))):
        if chrom not in scores:
            scores[chrom] = BinnedArray()
        scores[chrom][pos] = val

        # Status
        if i % 10000 == 0:
            print(i, "scores processed")

    for chr in scores.keys():
        out = open(chr, "w")
        scores[chr].to_file(out)
        out.close()


if __name__ == "__main__":
    main()
