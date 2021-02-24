#!/usr/bin/env python

"""
Read a maf from stdin and print the chromosome number for each alignment. It
searches for 'chr' in each alignment block src, and may not be robust if other
src formats are used.

NOTE: See 'align_print_template.py' for a more general variation of this
      program.

usage: %prog refindex [options]
"""

import sys

from bx.align import maf
from bx.cookbook import doc_optparse


def __main__():

    # Parse command line arguments
    options, args = doc_optparse.parse(__doc__)

    try:
        refindex = int(args[0])
    except Exception:
        doc_optparse.exit()

    maf_reader = maf.Reader(sys.stdin)

    for m in maf_reader:
        c = m.components[refindex].src
        print(c[c.rfind("chr") + 3:])


if __name__ == "__main__":
    __main__()
