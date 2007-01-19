#!/usr/bin/env python

"""
Read a file and print each line for which 'expr' is true. Expressions
can use variables "line" (the current line) and "fields" (line split
on '\t').

TODO: Does table_filter obviate this?

usage: %prog "expr" < file
"""

import psyco_full
import sys

def __main__():

    # Parse command line arguments

    expr = sys.argv[1]
    # Compile expression for SPEED
    expr = compile( expr, '<expr arg>', 'eval' )

    for line in sys.stdin:
        if bool( eval( expr, { "line": line, "fields": line.split("\t") } ) ):
            print line,

if __name__ == "__main__": __main__()
