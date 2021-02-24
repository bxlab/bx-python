#!/usr/bin/env python

"""
Read a MAF from standard input and print the score of each block. It can
optionally recalculate each score using the hox70 matrix, and normalize the
score by the number of columns in the alignment.

TODO: Should be able to read an arbitrary scoring matrix.

usage: %prog [options]
   -r, --recalculate: don't use the score from the maf, recalculate (using hox70 matrix)
   -l, --lnorm: divide (normalize) score by alignment text length
"""

import sys

from bx.align import maf, score
from bx.cookbook import doc_optparse


def main():

    # Parse command line arguments
    options, args = doc_optparse.parse(__doc__)

    try:
        lnorm = bool(options.lnorm)
        recalculate = bool(options.recalculate)
    except Exception:
        doc_optparse.exit()

    hox70 = score.build_scoring_scheme("""  A    C    G    T
                                      91 -114  -31 -123
                                    -114  100 -125  -31
                                     -31 -125  100 -114
                                    -123  -31 -114   91 """, 400, 30, default=0)

    maf_reader = maf.Reader(sys.stdin)

    for m in maf_reader:
        if m.text_size == 0:
            print("NA")
            continue
        s = m.score
        # Recalculate?
        if recalculate:
            s = hox70.score_alignment(m)
        # Normalize?
        if lnorm:
            s = s / m.text_size
        # Print
        print(s)


if __name__ == "__main__":
    main()
