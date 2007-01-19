#!/usr/bin/env python2.3

"""
Reads a list of intervals from `interval_file`, build an intersector and pickle it
to `pickle_file`.  

usage: %prog interval_file pickle_file
"""

import psyco_full

import from bx.cookbook import doc_optparse

from bx import intervals
import pickle
import sys


def __main__():

    # Parse Command Line

    options, args = doc_optparse.parse( __doc__ )

    try:
        range_filename = args[ 0 ]
        pickle_filename = args[ 1 ]
    except:
        doc_optparse.exit()

    # Load Intervals

    intersecter = intervals.Intersecter()
    for line in file( range_filename ):
        fields = line.split()
        intersecter.add_interval( intervals.Interval( int( fields[0] ), int( fields[1] ) ) )

    # Pickle
    
    out = open( pickle_filename, 'w' )
    pickle.dump( intersecter, out, pickle.HIGHEST_PROTOCOL )
    out.close()
    
if __name__ == "__main__": __main__()
