#!/usr/bin/env python

"""
Read some whitespace separated data from stdin and pretty print it so that
the columns line up.
"""

import sys


def main():
    pad = "\t"
    align = None
    if len(sys.argv) > 1:
        pad = " " * int(sys.argv[1])
    if len(sys.argv) > 2:
        align = sys.argv[2]
    rows = [line.split() for line in sys.stdin]
    print_tabular(rows, pad, align)


def print_tabular(rows, pad, align=None):
    if len(rows) == 0:
        return ""
    lengths = [len(col) for col in rows[0]]
    for row in rows[1:]:
        for i in range(0, len(row)):
            lengths[i] = max(lengths[i], len(row[i]))
    rval = ""
    for row in rows:
        rval = ""
        for i in range(0, len(row)):
            if align and align[i] == "l":
                rval += row[i].ljust(lengths[i])
            else:
                rval += row[i].rjust(lengths[i])
            rval += pad
        print(rval)


main()
