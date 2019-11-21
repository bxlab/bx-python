#!/usr/bin/env python

"""
Reads a list of intervals and a maf. Produces a new maf containing the
blocks or parts of blocks in the original that overlapped the intervals.

It is assumed that each file `maf_fname` has a corresponding `maf_fname`.index
file.

NOTE: If two intervals overlap the same block it will be written twice. With
      non-overlapping intervals and --chop this is never a problem.

NOTE: Intervals are origin-zero, half-open.  For example, the interval 100,150
      is 50 bases long, and there are 100 bases to its left in the sequence.

NOTE: Intervals are relative to the + strand, regardless of the strands in
      the alignments.


WARNING: bz2/bz2t support and file cache support are new and not as well
         tested.

usage: %prog maf_fname1 maf_fname2 ... [options] < interval_file
   -m, --mincols=0: Minimum length (columns) required for alignment to be output
   -c, --chop:       Should blocks be chopped to only portion overlapping (no by default)
   -s, --src=s:      Use this src for all intervals
   -p, --prefix=p:   Prepend this to each src before lookup
   -d, --dir=d:      Write each interval as a separate file in this directory
   -S, --strand:     Strand is included as an additional column, and the blocks are reverse complemented (if necessary) so that they are always on that strand w/r/t the src species.
   -C, --usecache:   Use a cache that keeps blocks of the MAF files in memory (requires ~20MB per MAF)
"""

import os
import sys

import bx.align.maf
from bx.cookbook import doc_optparse


def main():
    # Parse Command Line
    options, args = doc_optparse.parse(__doc__)
    try:
        maf_files = args
        if options.mincols:
            mincols = int(options.mincols)
        else:
            mincols = 0
        if options.src:
            fixed_src = options.src
        else:
            fixed_src = None
        if options.prefix:
            prefix = options.prefix
        else:
            prefix = None
        if options.dir:
            dir = options.dir
        else:
            dir = None
        chop = bool(options.chop)
        do_strand = bool(options.strand)
        use_cache = bool(options.usecache)
    except Exception:
        doc_optparse.exit()
    # Open indexed access to mafs
    index = bx.align.maf.MultiIndexed(maf_files, keep_open=True, parse_e_rows=True, use_cache=use_cache)
    # Start MAF on stdout
    if dir is None:
        out = bx.align.maf.Writer(sys.stdout)
    # Iterate over input ranges
    for line in sys.stdin:
        strand = None
        fields = line.split()
        if fixed_src:
            src, start, end = fixed_src, int(fields[0]), int(fields[1])
            if do_strand:
                strand = fields[2]
        else:
            src, start, end = fields[0], int(fields[1]), int(fields[2])
            if do_strand:
                strand = fields[3]
        if prefix:
            src = prefix + src
        # Find overlap with reference component
        blocks = index.get(src, start, end)
        # Open file if needed
        if dir:
            out = bx.align.maf.Writer(open(os.path.join(dir, "%s:%09d-%09d.maf" % (src, start, end)), 'w'))
        # Write each intersecting block
        if chop:
            for block in blocks:
                for ref in block.get_components_by_src(src):
                    slice_start = max(start, ref.get_forward_strand_start())
                    slice_end = min(end, ref.get_forward_strand_end())
                    if slice_end <= slice_start:
                        continue
                    sliced = block.slice_by_component(ref, slice_start, slice_end)
                    # If the block is shorter than the minimum allowed size, stop
                    if mincols and (sliced.text_size < mincols):
                        continue
                    # If the reference component is empty, don't write the block
                    if sliced.get_component_by_src(src).size < 1:
                        continue
                    # Keep only components that are not empty
                    sliced.components = [c for c in sliced.components if c.size > 0]
                    # Reverse complement if needed
                    if strand is not None and ref.strand != strand:
                        sliced = sliced.reverse_complement()
                    # Write the block
                    out.write(sliced)
        else:
            for block in blocks:
                out.write(block)
        if dir:
            out.close()
    # Close output MAF
    out.close()
    index.close()


if __name__ == "__main__":
    main()
