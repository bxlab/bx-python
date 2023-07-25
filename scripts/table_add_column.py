#!/usr/bin/env python

"""
Tool for adding a column to a table. Expressions for the column are similar
to those supported by table_filter.py

usage: %prog expression colname < table
    -H, --header:    keep header in output
    -C, --comments:  keep comments in output
"""

import sys

import bx.tabular.io
from bx.cookbook import doc_optparse


def __main__():
    # Parse command line arguments
    options, args = doc_optparse.parse(__doc__)
    try:
        keep_header = bool(options.header)
        keep_comments = bool(options.comments)
        expr = args[0]
        colname = args[1]
    except Exception:
        doc_optparse.exception()

    # Compile expression for SPEED
    if expr:
        expr = compile(expr, "<expr arg>", "eval")

    for element in bx.tabular.io.Reader(sys.stdin):
        if isinstance(element, bx.tabular.io.Header):
            if keep_header:
                print(str(element) + "\t" + colname)
        elif isinstance(element, bx.tabular.io.Comment):
            if keep_comments:
                print(element)
        else:
            val = eval(expr, dict(row=element))
            print(str(element) + "\t" + str(val))


if __name__ == "__main__":
    __main__()
