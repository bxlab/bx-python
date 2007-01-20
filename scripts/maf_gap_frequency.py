#!/usr/bin/env python

"""
Read a MAF from standard input and print the fraction of gap columns in 
each block.

usage: %prog < maf > out
"""

from __future__ import division

import sys
import bx.align.maf


def main():    
    for m in bx.align.maf.Reader( sys.stdin ):  
        gaps = 0        
        for col in m.column_iter():
            if ( '-' in col ): gaps += 1          
        print gaps / m.text_size

if __name__ == "__main__": main()
