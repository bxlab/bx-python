#!/usr/bin/env python

"""
Read data in UCSC wiggle format and write it to an "array tree" file.

usage: %prog array_length output.array_tree < input.wig
"""


import sys

from bx.arrays.array_tree import (
    array_tree_dict_from_reader,
    FileArrayTreeDict,
)
from bx.arrays.wiggle import WiggleReader


def main():
    sizes_fname = sys.argv[1]
    out_fname = sys.argv[2]

    sizes = {}
    for line in open(sizes_fname):
        fields = line.split()
        sizes[fields[0]] = int(fields[1])

    # Fill array from wiggle
    d = array_tree_dict_from_reader(WiggleReader(sys.stdin), sizes)

    for value in d.values():
        value.root.build_summary()

    with open(out_fname, "w") as f:
        FileArrayTreeDict.dict_to_file(d, f)


if __name__ == "__main__":
    main()
