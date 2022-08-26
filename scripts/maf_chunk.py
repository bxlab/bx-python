#!/usr/bin/env python

"""
Read a MAF from stdin and break into several new mafs containing no more than
`chunk_size` columns. The new mafs will be written to `out_dir` along with a
file "intervals.txt" specifying the range covered by each new maf file. A
probability for writing each chunk can optionally be specified, resulting in
a random fraction of chunks from the input MAF being produced.

usage: %prog [options] chunk_size out_dir < maf
  --prob: probability of writing versus skipping each chunk.
"""

import random
import sys
from optparse import OptionParser

import numpy as np

import bx.align.maf

INF = np.inf


def __main__():
    parser = OptionParser("usage: %prog chunk_size out_dir")
    parser.add_option("--prob", action="store", default=None, type="float", help="Probability of writing a given chunk")

    (options, args) = parser.parse_args()

    chunk_size = int(args[0])
    out_dir = args[1]
    prob = options.prob

    maf_reader = bx.align.maf.Reader(sys.stdin, parse_e_rows=True)

    maf_writer = None

    count = 0
    current_chunk = -1

    chunk_min = INF
    chunk_max = 0

    write_current_chunk = True

    interval_file = open("%s/intervals.txt" % out_dir, "w")

    for m in maf_reader:
        if not maf_writer or count + m.text_size > chunk_size:
            current_chunk += 1
            # Finish the last chunk
            if maf_writer:
                maf_writer.close()
                interval_file.write(f"{chunk_min} {chunk_max}\n")
                chunk_min = INF
                chunk_max = 0
            # Decide if the new chunk will be written
            if prob:
                write_current_chunk = bool(random.random() <= prob)
            else:
                write_current_chunk = True
            if write_current_chunk:
                maf_writer = bx.align.maf.Writer(open("%s/%09d.maf" % (out_dir, current_chunk), "w"))
            else:
                maf_writer = None
            count = 0
        if maf_writer:
            maf_writer.write(m)
        # count += m.text_size
        count += m.components[0].size
        chunk_min = min(chunk_min, m.components[0].start)
        chunk_max = max(chunk_max, m.components[0].end)

    if maf_writer:
        maf_writer.close()
        interval_file.write(f"{chunk_min} {chunk_max}\n")

    interval_file.close()


if __name__ == "__main__":
    __main__()
