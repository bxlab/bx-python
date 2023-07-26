#!/usr/bin/env python

"""
Read scores in "wiggle" format from `score_file` and intervals in "bed" format
from `interval_file` and print all scores overlapping intervals.

TODO: Support binned array format scores also.

usage: %prog score_file interval_file [out_file]
"""

import sys

import bx.wiggle
from bx import misc
from bx.binned_array import BinnedArray
from bx.cookbook import doc_optparse


def read_scores(f):
    scores_by_chrom = dict()
    for chrom, pos, val in bx.wiggle.Reader(f):
        if chrom not in scores_by_chrom:
            scores_by_chrom[chrom] = BinnedArray()
        scores_by_chrom[chrom][pos] = val
    return scores_by_chrom


def main():
    # Parse command line
    options, args = doc_optparse.parse(__doc__)
    try:
        score_file = open(args[0])
        interval_file = open(args[1])
        if len(args) > 2:
            out_file = open(args[2], "w")
        else:
            out_file = sys.stdout
    except Exception:
        doc_optparse.exit()

    scores_by_chrom = read_scores(misc.open_compressed(sys.argv[1]))
    for line in open(sys.argv[2]):
        fields = line.split()
        chrom, start, stop = fields[0], int(fields[1]), int(fields[2])
        if chrom in scores_by_chrom:
            ba = scores_by_chrom[chrom]
            scores = [ba[i] for i in range(start, stop)]
        else:
            scores = []
        print(" ".join(fields), " ".join(map(str, scores)), file=out_file)

    score_file.close()
    interval_file.close()
    out_file.close()


if __name__ == "__main__":
    main()
