"""
Complement a set of intervals.
"""

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *
from bx.bitset import MAX

def complement( reader, lens ):
    # Handle any ValueError, IndexError and OverflowError exceptions that may be thrown when
    # the bitsets are being created by skipping the problem lines
    complement_reader = BitsetSafeReaderWrapper( reader, lens=lens )
    bitsets = complement_reader.binned_bitsets( upstream_pad=0, downstream_pad=0, lens=lens )
    # NOT them all
    for key, value in bitsets.items():
        value.invert()
    # Read remaining intervals and subtract
    for chrom in bitsets:
        bitset = bitsets[chrom]
        out_intervals = bits_set_in_range( bitset, 0, lens.get( chrom, MAX ) )
        try:
            # Write the intervals
            for start, end in out_intervals:
                fields = ["."  for x in range(max(complement_reader.chrom_col, complement_reader.start_col, complement_reader.end_col)+1)]
                # default the column to a + if it exists
                if complement_reader.strand_col < len( fields ) and complement_reader.strand_col >= 0:
                    fields[complement_reader.strand_col] = "+"
                fields[complement_reader.chrom_col] = chrom
                fields[complement_reader.start_col] = start
                fields[complement_reader.end_col] = end
                new_interval = GenomicInterval(complement_reader, fields, complement_reader.chrom_col, complement_reader.start_col, complement_reader.end_col, complement_reader.strand_col, "+")
                yield new_interval
        except IndexError, e:
            complement_reader.skipped += 1
            # no reason to stuff an entire bad file into memmory
            if complement_reader.skipped < 10:
                complement_reader.skipped_lines.append( ( complement_reader.linenum, complement_reader.current_line, str( e ) ) )
            continue


# def main():
#     # test it all out
#     f1 = fileinput.FileInput("dataset_7.dat")
#     g1 = GenomicIntervalReader(f1)
#     for interval in complement(g1,{"chr":16000000}):
#         print "\t".join(interval)
# 
# if __name__ == "__main__":
#     main()
