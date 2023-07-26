#!/usr/bin/env python

"""
Randomly shuffle the columns of each block of a maf file. Note that this does
not change any other features of the maf block, thus the text of each row no
longer will match the sequence refered to by the other row attributes!

usage: %prog < maf > maf
"""

import sys

from bx import align
from bx.align import maf


def __main__():
    maf_reader = maf.Reader(sys.stdin, parse_e_rows=True)
    maf_writer = maf.Writer(sys.stdout)

    for m in maf_reader:
        align.shuffle_columns(m)

        maf_writer.write(m)


if __name__ == "__main__":
    __main__()
