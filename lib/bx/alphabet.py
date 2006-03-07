from numarray import *

import seq_numarray

DNA_BASE=5

## FIXME: Is DNA_BASE really neccesary? This should be able to map from
##        any integer alphabet to any other.

class Mapping( object ):
    def __init__( self, f=None ):
        if f: self.read_from_file( f )

    def read_from_file( self, f ):
        align_count = None
        max_symbol = 0
        reverse_map = {}
        for line in f:
            ( key, value ) = line.split()
            if not align_count: 
                align_count = len( key )
                self.table = zeros( DNA_BASE ** align_count )
            else:
                assert align_count == len( key )
            index = seq_numarray.DNA.translate_alignment_column( key )
            self.table[ index ] = int( value )
            reverse_map.setdefault( int( value ), [] ).append( key )
            max_symbol = max( self.table[ index ], max_symbol )
        self.align_count = align_count
        self.symbol_count = max_symbol + 1
        self.reverse_table = [ reverse_map[ index ] for index in range( 0, self.symbol_count ) ]

    def collapse( self, a, b ):
        copy = Mapping()
        copy.align_count = self.align_count
        copy.symbol_count = self.symbol_count - 1
        table = self.table.copy()
        for i in range( len( table ) ):
            if table[i] == b: table[i] = a
            elif table[i] == copy.symbol_count: table[i] = b
        copy.table = table
        return copy
        
    def translate_alignment( self, seqs ):
        return self.translate( seq_numarray.DNA.translate_alignment( seqs ) )

    def translate( self, seq ):
        return take( self.table, seq )

    def reverse( self, seq ):
        return [ self.reverse_table[ index ] for index in seq ]
