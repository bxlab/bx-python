#!/usr/bin/env python

"""
Mask out potential CpG sites from a maf. Restricted or inclusive definition
of CpG sites can be used. The total fraction masked is printed to stderr.

usage: %prog < input > output
    -m, --mask=N: Character to use as mask ('?' is default)
    -r, --restricted: Use restricted definition of CpGs
"""

import sys

import bx.align
import bx.align.maf
import bx.align.sitemask.cpg
from bx.cookbook import doc_optparse


def main():
    options, args = doc_optparse.parse(__doc__)
    try:
        if options.mask:
            mask = options.mask
        else:
            mask = "?"
    except Exception:
        doc_optparse.exception()

    reader = bx.align.maf.Reader(sys.stdin)
    writer = bx.align.maf.Writer(sys.stdout)

    if options.restricted:
        cpgfilter = bx.align.sitemask.cpg.Restricted(mask=mask)
    else:
        cpgfilter = bx.align.sitemask.cpg.Inclusive(mask=mask)
    cpgfilter.run(reader, writer.write)

    print(str(float(cpgfilter.masked) / float(cpgfilter.total) * 100) + "% bases masked.", file=sys.stderr)


if __name__ == "__main__":
    main()
