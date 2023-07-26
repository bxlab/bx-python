#!/usr/bin/env python

"""
Filter maf blocks for presence of wildcard columns. Blocks must meet the
criteria of having at least `min_good` columns, each of which has more than
`min_species` rows that are NOT wildcard bases ('*').

TODO: Allow specifying the character of the wildcard base.

usage: %prog min_good min_species < maf > maf
"""

import sys

from bx.align import maf


def main():
    min_good = int(sys.argv[1])
    min_species = int(sys.argv[2])

    maf_reader = maf.Reader(sys.stdin, parse_e_rows=True)
    maf_writer = maf.Writer(sys.stdout)

    for m in maf_reader:
        good = 0
        for col in m.column_iter():
            if col.count("*") <= min_species:
                good += 1
        if good >= min_good:
            maf_writer.write(m)


if __name__ == "__main__":
    main()
