#!/usr/bin/env python2.3

from numarray import *
from numarray.ieeespecial import *
from cookbook.attribute import *

import string

class SeqTranslater( object ):
    attribute( table=None, base=0 )
    """Encapsulates a mapping from chars to integers"""
    def __init__( self, mapping ):
        """Initialize translation table from a list in which the 
           characters of each list item map to the same number"""
        table = zeros( 256 ) - 1
        reverse_mapping = []
        for value, chars in enumerate( mapping ):
            for ch in chars:
                put( table, array( ch, 'b' ), value )
            reverse_mapping.append( chars[ 0 ] )
        self.table = table 
        self.reverse_mapping = reverse_mapping
        self.base = len( mapping )
    def translate( self, seq ):
        """Convert a character sequence to a single integer array"""
        int_seq = take( self.table, array( seq, 'b' ) )
        if -1 in int_seq: raise "Invalid character in sequence"
        return int_seq
    def translate_alignment( self, seqs ):   
        """Convert the rows of a multiple alignment to a single integer array"""
        if len( seqs ) < 1: return None
        rval = zeros( len( seqs[ 0 ] ) )
        factor = 1
        for seq in seqs:
            seq_ints = self.translate( seq )
            rval += ( seq_ints * factor )
            factor *= self.base
        return rval
    def translate_alignment_column( self, col, allow_invalid=False ):
        value = 0
        factor = 1
        for ch in col:
            row_value = self.table[ int( array( ch, 'b' )[0] ) ]
            if row_value == -1: 
                if allow_invalid: return -1
                else: raise "Invalid character"
            value += ( row_value * factor )
            factor *= self.base
        return value
    def reverse_alignment_column( self, align_count, value ):
        col = []
        for i in range( 0, align_count ):
            index = ( value / self.base ** i ) % self.base 
            col.append( self.reverse_mapping[ index ] )
        return string.join( col, '' )

DNA = SeqTranslater( ( 'Aa', 'Cc', 'Gg', 'Tt', '-NnBb' ) )
