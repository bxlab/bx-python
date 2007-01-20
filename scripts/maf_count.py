#!/usr/bin/env python2.3

"""
Read a MAF from standard input and print counts of alignments, bases, or 
columns. 

usage: %prog [options]
   -c, --cols: count alignment columns rather than number of alignments
   -b, --bases: count bases in first species rather than number of alignments
   -s, --skip=N: when counting bases, skip this base
   -e, --each: print a count for each alignment rather than whole file
   -r, --ref=N: reference sequence (first by default, 0..n)
"""

from bx.cookbook import doc_optparse
import sys

import bx.align.maf

def __main__():

    options, args = doc_optparse.parse( __doc__ )

    try:
        if options.cols: action = "cols"
        elif options.bases: action = "bases"
        else: action = "aligns"
        print_each = bool( options.each )
        if options.ref: ref = int( options.ref )
        else: ref = 0
        if options.skip: skip = options.skip
        else: skip = None
    except:
        doc_optparse.exit()

    maf_reader = bx.align.maf.Reader( sys.stdin )

    count = 0

    for m in maf_reader:
        
        if action == "aligns": 
            count += 1
        elif action == "cols": 
            count += m.text_size
        elif action == "bases":
            if skip:
                count += ( m.components[ref].size - m.components[ref].text.count( skip ) )
            else:
                count += m.components[ref].size

        if print_each: 
            print count
            count = 0

    if not print_each: print count

if __name__ == "__main__": __main__()
