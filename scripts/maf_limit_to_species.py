#!/usr/bin/env python

"""
Read a maf file from stdin and write out a new maf with only blocks having all
of the required in species, after dropping any other species and removing
columns containing only gaps.

usage: %prog species,species2,... < maf
"""

import sys

import bx.align.maf


def main():
    species = sys.argv[1].split(",")

    maf_reader = bx.align.maf.Reader(sys.stdin, parse_e_rows=True)
    maf_writer = bx.align.maf.Writer(sys.stdout)

    for m in maf_reader:
        new_components = []
        for comp in m.components:
            if comp.src.split(".")[0] in species:
                new_components.append(comp)
        m.components = new_components
        m.remove_all_gap_columns()
        if len(m.components) > 1:
            maf_writer.write(m)

    maf_reader.close()
    maf_writer.close()


if __name__ == "__main__":
    main()
