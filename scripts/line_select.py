#!/usr/bin/env python2.3

"""
Read a feature file containing a 0 or 1 on each line, output 
all lines from stdin for which the feature file is 1

usage: %prog feature_file < ...
"""

import psyco_full

import sys

def __main__():

    feature_file = sys.argv[1]

    if len( sys.argv ) > 2:
        match = int( sys.argv[2] )
    else:
        match = 1
    
    feature_vector = [ int( line ) for line in file( feature_file ) ]

    for index, line in enumerate( sys.stdin ):
        if feature_vector[ index ] == match: print line,

if __name__ == "__main__": __main__()
