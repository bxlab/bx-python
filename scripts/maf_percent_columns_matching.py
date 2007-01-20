#!/usr/bin/env python2.3

"""
Read a PAIRWISE maf from stdin and print the fraction of columns whose bases 
match for each alignment.

TODO: generalize for more than two speceis.

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
            if a == b:
                match += 1
            total += 1

        print match / total


if __name__ == "__main__": __main__()
