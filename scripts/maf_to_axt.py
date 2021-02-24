#!/usr/bin/env python
"""
Application to convert MAF file to AXT file, projecting to any two species.
Reads a MAF file from standard input and writes an AXT file to standard out;
some statistics are written to standard error.  The user must specify the
two species of interest.

usage: %prog primary_species secondary_species < maf_file > axt_file
"""

__author__ = "Bob Harris (rsharris@bx.psu.edu)"

import copy
import sys

import bx.align.axt
import bx.align.maf


def usage(s=None):
    message = """
maf_to_axt primary_species secondary_species < maf_file > axt_file
"""
    if s is None:
        sys.exit(message)
    else:
        sys.exit(f"{s}\n{message}")


def main():
    primary = None
    secondary = None

    args = sys.argv[1:]
    while len(args) > 0:
        arg = args.pop(0)
        val = None
        fields = arg.split("=", 1)
        if len(fields) == 2:
            arg = fields[0]
            val = fields[1]
            if val == "":
                usage("missing a value in %s=" % arg)

        if primary is None and val is None:
            primary = arg
        elif secondary is None and val is None:
            secondary = arg
        else:
            usage("unknown argument: %s" % arg)

    if primary is None:
        usage("missing primary species")

    if secondary is None:
        usage("missing secondary species")

    # read the alignments and other info

    out = bx.align.axt.Writer(sys.stdout)

    axtsRead = 0
    mafsWritten = 0
    for mafBlock in bx.align.maf.Reader(sys.stdin):
        axtsRead += 1

        p = mafBlock.get_component_by_src_start(primary)
        if p is None:
            continue
        s = mafBlock.get_component_by_src_start(secondary)
        if s is None:
            continue

        axtBlock = bx.align.Alignment(mafBlock.score, mafBlock.attributes)
        axtBlock.add_component(clone_component(p))
        axtBlock.add_component(clone_component(s))

        remove_mutual_gaps(axtBlock)
        if axtBlock.text_size == 0:
            continue

        out.write(axtBlock)
        mafsWritten += 1

    sys.stderr.write("%d blocks read, %d written\n" % (axtsRead, mafsWritten))


def clone_component(c):
    return bx.align.Component(c.src, c.start, c.size, c.strand, c.src_size, copy.copy(c.text))


def remove_mutual_gaps(block):

    if len(block.components) == 0:
        return

    nonGaps = []

    for c in block.components:
        for ix in range(0, block.text_size):
            if ix not in nonGaps and c.text[ix] != "-":
                nonGaps.append(ix)

    nonGaps.sort()

    for c in block.components:
        c.text = "".join([c.text[ix] for ix in nonGaps])

    block.text_size = len(nonGaps)


if __name__ == "__main__":
    main()
