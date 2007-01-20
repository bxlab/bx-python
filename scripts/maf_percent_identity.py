#!/usr/bin/env python2.3

"""
Read a PAIRWISE maf from stdin and print the percent identity of each
alignment, where percent identity is defined as the number of matching columns
over the number of aligned (non-gap) columns.

TODO: Generalize for more than two species

usage: %prog < maf > out
"""

from __future__ import division

import sys

import psyco_full

from bx.align import maf


def __main__():

    maf_reader = maf.Reader( sys.stdin )

    for m in maf_reader:
        match = 0
        total = 0
        for i in range( 0, m.text_size ):
            a = m.components[0].text[i].lower() 
            b = m.components[1].text[i].lower()            
            if a == '-' or b == '-': 
                continue
            elif a == b:
                match += 1
            total += 1

        print match / total


if __name__ == "__main__": __main__()
