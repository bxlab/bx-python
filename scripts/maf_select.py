#!/usr/bin/env python

"""
Read a feature file containing a 0 or 1 on each line, output
all mafs whose index in maf_file corresponds to a row having a 1

usage: %prog feature_file < maf_file
"""

import sys

import bx.align.maf


def __main__():

    feature_file = sys.argv[1]

    if len(sys.argv) > 2:
        match = int(sys.argv[2])
    else:
        match = 1

    feature_vector = [int(line) for line in open(feature_file)]

    maf_reader = bx.align.maf.Reader(sys.stdin)
    maf_writer = bx.align.maf.Writer(sys.stdout)

    index = 0

    for m in maf_reader:
        if feature_vector[index] == match:
            maf_writer.write(m)
        index += 1


if __name__ == "__main__":
    __main__()
