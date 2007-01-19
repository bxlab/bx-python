#!/usr/bin/env python

import psyco; psyco.profile()

from numarray import *
from numarray.ieeespecial import *

from tables import *

from bx.misc import *

from bx.cookbook.progress_bar import ProgressBar

def main():

    lengths = read_lengths( open( sys.argv[1] ) )
    out_filename = sys.argv[2]
    
    h5out = openFile( out_filename, mode = "a", filters=Filters( complevel=7 ) )

    for fname in sys.argv[3:]:
        chr = fname.split( '/' )[-1].split( '.' )[0]
        # Initialize array to NaN
        length = lengths[ chr ]
        a = zeros( length, type=Float32 )
        a[:] = nan
        bar = ProgressBar( 0, length, 80 )
        # Read scores into array
        print >> sys.stderr, 'Reading scores from ' + fname
        for i, line in enumerate( open_compressed( fname ) ):
            fields = line.split()
            pos = int( fields[ 0 ] ) 
            a[ pos ] = float( fields[1] )
            if ( i % 1000 ) == 0: bar.update_and_print( pos )

        # Write to h5 file
        print >> sys.stderr, 'Adding result to hdf file as /' + chr
        atom = Float32Atom( shape=(0,) )
        arr = h5out.createEArray( h5out.root, chr, atom, "Scores for " + chr )
        arr.append( a )
    
    h5out.close()

def read_lengths( f ):
    lengths = {}
    for line in f:
        fields = line.split()
        lengths[ fields[0].split('.')[0] ] = int( fields[1] )
    return lengths    
        
main()
