"""
Rudimentary support for PHAST's tree model file format

ALPHABET: A C G T - 
ORDER: 0
SUBST_MOD: HKY85+Gap
TRAINING_LNL: -8760441383.526642
BACKGROUND: 0.265517 0.184702 0.184801 0.265947 0.099033 
RATE_MAT:
...
"""

from scipy import *
from scipy.linalg import *

from cookbook import *

class TreeModel:
    def __init__( self ):
        self.alphabet = None
        self.radix = 0
        self.background = None
        self.tree = None
        self.matrix = None
    def matrix_for_time( self, t ):
        return expm( self.matrix * t )
    matrix_for_time = cachedmethod( matrix_for_time )

def from_file( f ):
    input = iter( f )
    tm = TreeModel()
    for line in input:
        if line.startswith( "ALPHABET:" ):
            tm.alphabet = tuple( line.split()[1:] )    
            tm.radix = len( tm.alphabet )
        if line.startswith( "BACKGROUND:" ):
            tm.background = tuple( map( float, line.split()[1:] ) )
        if line.startswith( "TREE:" ):
            tm.tree = line[6:].strip() 
        if line.startswith( "RATE_MAT:" ):
            matrix = zeros( (tm.radix,tm.radix), Float )
            for i in range( len( tm.alphabet ) ):
                matrix[i] = map( float, input.next().split() )
            tm.matrix = matrix
    return tm

if __name__ == "__main__":
    import sys
    tm = from_file( open( sys.argv[1] ) )
    print tm
    print tm.alphabet
    print tm.background
    print tm.matrix
    print tm.matrix_for_time( 1 )
