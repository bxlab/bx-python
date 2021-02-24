#!/usr/bin/env python

"""
Tool for filtering a tabular data file. Fields are separated by tabs, the
header line is denoted by a '#' in the first byte, comments are denoted by
a '#' at the start of any subsequent line.

Expressions can use column names as well as numbers. The -c options allows
cutting, again using field name or numbers.

usage: %prog expression < table
    -H, --header:       keep header in output
    -C, --comments:     keep comments in output
    --force-header:     assume the first line is a header even if it does not start with "#"
    -c, --cols=1,2:     names or indexes of columns to keep
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
        cols = []
        if options.cols:
            for c in options.cols.split(','):
                try:
                    v = int(c)
                except ValueError:
                    v = c
                cols.append(v)
        if len(args) > 0:
            expr = args[0]
        else:
            expr = None
        if options.force_header:
            force_header = bx.tabular.io.FIRST_LINE_IS_HEADER
        else:
            force_header = None
    except Exception:
        doc_optparse.exception()

    # Compile expression for SPEED
    if expr:
        expr = compile(expr, '<expr arg>', 'eval')

    for element in bx.tabular.io.TableReader(sys.stdin, force_header=force_header):
        if isinstance(element, bx.tabular.io.Header):
            if keep_header:
                if cols:
                    print("#" + "\t".join(element[c] for c in cols))
                else:
                    print(element)
        elif isinstance(element, bx.tabular.io.Comment):
            if keep_comments:
                print(element)
        else:
            if expr is None or bool(eval(expr, dict(row=element))):
                if cols:
                    print("\t".join(element[c] for c in cols))
                else:
                    print(element)


if __name__ == "__main__":
    __main__()
