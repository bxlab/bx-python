#!/usr/bin/env python

"""
Read a maf file from stdin and write out a new maf with only blocks having all of
the passed in species, after dropping any other species and removing columns
containing only gaps. By default this will attempt to fuse together any blocks
which are adjacent after the unwanted species have been dropped.

usage: %prog species1 species2 ... < maf
    -n, --nofuse: Don't attempt to join blocks, just remove rows.
"""

import sys

import bx.align.maf
from bx.align.tools.fuse import FusingAlignmentWriter
from bx.align.tools.thread import get_components_for_species
from bx.cookbook import doc_optparse


def main():
    options, args = doc_optparse.parse(__doc__)

    try:
        species = args
        # Allow a comma separated list, TODO: allow a newick format tree
        if len(species) == 1 and "," in species[0]:
            species = species[0].split(",")
        fuse = not (bool(options.nofuse))
    except Exception:
        doc_optparse.exit()

    maf_reader = bx.align.maf.Reader(sys.stdin, parse_e_rows=True)
    maf_writer = bx.align.maf.Writer(sys.stdout)

    if fuse:
        maf_writer = FusingAlignmentWriter(maf_writer)

    for m in maf_reader:
        new_components = get_components_for_species(m, species)
        if new_components:
            m.components = new_components
            m.score = 0.0
            m.remove_all_gap_columns()
            maf_writer.write(m)

    maf_reader.close()
    maf_writer.close()


if __name__ == "__main__":
    main()
