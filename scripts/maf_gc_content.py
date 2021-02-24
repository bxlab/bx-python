#!/usr/bin/env python

"""
Read a MAF from standard input and print average GC content of each alignment

usage: %prog < maf > out
"""

import sys

from bx.align import maf


def __main__():

    maf_reader = maf.Reader(sys.stdin)

    for m in maf_reader:
        gc = 0
        bases = 0
        for c in m.components:
            gc += c.text.count('G')
            gc += c.text.count('C')
            gc += c.text.count('g')
            gc += c.text.count('c')
            bases += (len(c.text) - c.text.count('-'))

        print(gc / bases)


if __name__ == "__main__":
    __main__()
