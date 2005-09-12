#!/usr/bin/env python

"""
usage: %prog score_file interval_file [out_file] [options] 
"""

from __future__ import division

import sys
import psyco_full
import bx.wiggle
from bx.binned_array import BinnedArray
from fpconst import isNaN
import cookbook.doc_optparse
import misc

def read_scores( f ):
    scores_by_chrom = dict()
    for chrom, pos, val in bx.wiggle.Reader( f ):
        if chrom not in scores_by_chrom:
            scores_by_chrom[chrom] = BinnedArray()
        scores_by_chrom[chrom][pos] = val
    return scores_by_chrom

def main():

    # Parse command line
    options, args = cookbook.doc_optparse.parse( __doc__ )
    try:
        score_file = open( args[0] )
        interval_file = open( args[1] )
        if len( args ) > 2:
            out_file = open( args[2], 'w' )
        else:
            out_file = sys.stdout
    except:
        cookbook.doc_optparse.exit()

    scores_by_chrom = read_scores( misc.open_compressed( sys.argv[1] ) )
    for line in open( sys.argv[2] ):
        fields = line.split()
        chrom, start, stop = fields[0], int( fields[1] ), int( fields[2] )
        total = 0
        count = 0
        min_score = 100000000
        max_score = -100000000
        for i in range( start, stop ):
            if chrom in scores_by_chrom and scores_by_chrom[chrom][i]:
                score = scores_by_chrom[chrom][i]
                if not isNaN( score ):
                    total += score
                    count += 1
                    max_score = max( score, max_score )
                    min_score = min( score, min_score )
        if count > 0:
            avg = total/count
        else:
            avg = "nan"
            min_score = "nan"
            max_score = "nan"
            
        print >> out_file, chrom, start, stop, avg, min_score, max_score

    score_file.close()
    interval_file.close()
    out_file.close()

if __name__ == "__main__": main()
