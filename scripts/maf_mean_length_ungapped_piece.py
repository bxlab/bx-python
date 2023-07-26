#!/usr/bin/env python

"""
Read a MAF from standard input and determine the mean length of ungapped pieces
in each block.

usage: %prog < maf > out
"""

import sys

import bx.align.maf


def main():
    for m in bx.align.maf.Reader(sys.stdin):
        ungapped_columns = 0
        ungapped_runs = 0
        in_ungapped = False

        for col in m.column_iter():
            is_gap = "-" in col
            if not is_gap:
                ungapped_columns += 1
            if in_ungapped and is_gap:
                ungapped_runs += 1
            in_ungapped = not is_gap
        if in_ungapped:
            ungapped_runs += 1

        print(ungapped_columns / ungapped_runs)


if __name__ == "__main__":
    main()
