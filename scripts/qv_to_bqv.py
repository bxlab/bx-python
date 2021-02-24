#!/usr/bin/env python

"""
Convert a qual (qv) file to several BinnedArray files for fast seek.
This script takes approximately 4 seconds per 1 million base pairs.

The input format is fasta style quality -- fasta headers followed by
whitespace separated integers.

usage: %prog qual_file output_file
"""

import fileinput
import sys

from bx.binned_array import BinnedArrayWriter


def main():
    args = sys.argv[1:]
    try:
        qual_file = args[0]
        output_file = args[1]
    except IndexError:
        print("usage: qual_file output_file")
        sys.exit()

    qual = fileinput.FileInput(qual_file)
    outfile = None
    outbin = None
    base_count = 0
    mega_count = 0
    region = ""

    for line in qual:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            # close old
            if outbin and outfile:
                print("\nFinished region " + region + " at " + str(base_count) + " base pairs.")
                outbin.finish()
                outfile.close()
            # start new file
            region = line.lstrip(">")
            outfname = output_file + "." + region + ".bqv"
            print("Writing region " + region + " to file " + outfname)
            outfile = open(outfname, "wb")
            outbin = BinnedArrayWriter(outfile, typecode='b', default=0)
            base_count = 0
            mega_count = 0
        else:
            if outfile and outbin:
                nums = line.split()
                for val in nums:
                    outval = int(val)
                    assert outval <= 255 and outval >= 0
                    outbin.write(outval)
                    base_count += 1
                if (mega_count * 1000000) <= base_count:
                    sys.stdout.write(str(mega_count)+" ")
                    sys.stdout.flush()
                    mega_count = base_count // 1000000 + 1
    if outbin and outfile:
        print("\nFinished region " + region + " at " + str(base_count) + " base pairs.")
        outbin.finish()
        outfile.close()


if __name__ == "__main__":
    main()
