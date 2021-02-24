#!/usr/bin/env python
"""
Create a bed file listing all the divergent sites between two specific species
in a maf.

usage: %prog maf_file reference_species_name other_species_name
"""

import sys

import bx.align.maf
import bx.bitset


def main():
    bitsets = {}
    maf = sys.argv[1]
    reference_sp, other_sp = sys.argv[2], sys.argv[3]

    for block in bx.align.maf.Reader(open(maf)):
        ref = block.get_component_by_src_start(reference_sp)
        other = block.get_component_by_src_start(other_sp)

        if not ref or not other:
            continue
        ref_chrom = ref.src.split('.')[1]
        ref_start = ref.start
        chrom_size = ref.get_src_size()

        if ref_chrom not in bitsets:
            bitsets[ref_chrom] = bx.bitset.BinnedBitSet(chrom_size)

        pos = ref_start
        for i, j in zip(ref.text.upper(), other.text.upper()):
            if i != '-':
                if i != j:  # mismatch
                    if i != 'N' and j != 'N' and j != '-':
                        # set if all valid chars
                        bitsets[ref_chrom].set(pos)
                pos += 1

    # bits --> bed file
    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while True:
            start = bits.next_set(end)
            if start == bits.size:
                break
            end = bits.next_clear(start)
            print("%s\t%d\t%d" % (chrom, start, end))


main()
