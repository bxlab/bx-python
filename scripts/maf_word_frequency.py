#!/usr/bin/env python

"""
Read a MAF and print counts and frequencies of all n-mers
(words composed on n consecutive alignment columns)

TODO: reconcile this and maf_mapping_word_frequency.py

usage: %prog n < maf_file
"""

import string
import sys

import psyco

from bx.align import maf

psyco.profile()


def __main__():

    motif_len = int(sys.argv[1])

    big_map = {}
    total = 0

    maf_reader = maf.Reader(sys.stdin)

    for m in maf_reader:
        texts = [c.text.upper() for c in m.components]
        for i in range(m.text_size - motif_len):
            motif = string.join([text[i: i + motif_len] for text in texts])
            if motif in big_map:
                big_map[motif] += 1
            else:
                big_map[motif] = 1
            total += 1

    items = sorted(zip(big_map.values(), big_map.keys()))
    items.reverse()

    for count, motif in items:
        print("%d\t%0.10f\t%s" % (count, count / total, motif))


if __name__ == "__main__":
    __main__()
