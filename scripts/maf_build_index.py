#!/usr/bin/env python

"""
Build an index file for a set of MAF alignment blocks.

If index_file is not provided maf_file.index is used.

usage: %prog maf_file index_file
    -s, --species=a,b,c: only index the position of the block in the listed species
"""

import os.path
from io import TextIOWrapper

import bx.align.maf
from bx import interval_index_file
from bx.cookbook import doc_optparse
from bx.misc.seekbzip2 import SeekableBzip2File
from bx.misc.seeklzop import SeekableLzopFile


def main():
    options, args = doc_optparse.parse(__doc__)

    try:
        maf_file = args[0]
        # If it appears to be a bz2 file, attempt to open with table
        if maf_file.endswith(".bz2"):
            table_file = maf_file + "t"
            if not os.path.exists(table_file):
                doc_optparse.exit("To index bz2 compressed files first "
                                  "create a bz2t file with bzip-table.")
            # Open with SeekableBzip2File so we have tell support
            maf_in = SeekableBzip2File(maf_file, table_file)
            # Strip .bz2 from the filename before adding ".index"
            maf_file = maf_file[:-4]
        elif maf_file.endswith(".lzo"):
            table_file = maf_file + "t"
            if not os.path.exists(table_file):
                doc_optparse.exit("To index lzo compressed files first "
                                  "create a lzot file with lzop_build_offset_table.")
            # Open with SeekableBzip2File so we have tell support
            maf_in = SeekableLzopFile(maf_file, table_file)
            # Strip .lzo from the filename before adding ".index"
            maf_file = maf_file[:-4]
        else:
            maf_in = open(maf_file, "rb")
        # Determine the name of the index file
        if len(args) > 1:
            index_file = args[1]
        else:
            index_file = maf_file + ".index"
        if options.species:
            species = options.species.split(",")
        else:
            species = None
    except Exception:
        doc_optparse.exception()

    maf_in = TextIOWrapper(maf_in, encoding="ascii")
    maf_reader = bx.align.maf.Reader(maf_in)

    indexes = interval_index_file.Indexes()

    # Need to be a bit tricky in our iteration here to get the 'tells' right
    while True:
        pos = maf_reader.file.tell()
        block = next(maf_reader)
        if block is None:
            break
        for c in block.components:
            if species is not None and c.src.split('.')[0] not in species:
                continue
            indexes.add(c.src, c.forward_strand_start, c.forward_strand_end, pos, max=c.src_size)

    out = open(index_file, 'wb')
    indexes.write(out)
    out.close()


if __name__ == "__main__":
    main()
