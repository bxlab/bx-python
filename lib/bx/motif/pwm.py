"""
Classes for working with position specific matrices.
"""

from numpy import *
from copy import copy

import _pwm

class BaseMatrix( object ):
    """
    Base class for position specific matrices. 
    """
    def __init__( self, alphabet=None, sorted_alphabet=None, 
                  char_to_index=None, values=None ):
        self.alphabet = alphabet
        self.sorted_alphabet = sorted_alphabet
        self.char_to_index = char_to_index
        self.values = values

    @classmethod
    def from_rows( Class, alphabet, rows ):
        """
        Create a new matrix for a sequence over alphabet `alphabet` taking 
        values from `rows` which is a list whose length is the width of the
        matrix, and whose elements are lists of values associated with each
        character (in the order those characters appear in alphabet). 
        """
        # Sorted alphabet
        sorted_alphabet = sorted( alphabet )
        # Character to index mapping (initialized to -1)
        char_to_index = zeros( (256), int16 ) - 1
        for i, ch  in enumerate( sorted_alphabet ):
            char_to_index[ ord(ch) ] = i
        # Array
        values = zeros( ( len( rows) , len( alphabet ) ), float32 )
        for i, row in enumerate( rows ):
            assert len( row ) == len( alphabet )
            for ch, val in zip( alphabet, row ):
                values[i, char_to_index[ord(ch)]] = val
        # Matrix
        matrix = Class()
        matrix.alphabet = alphabet
        matrix.sorted_alphabet = sorted_alphabet
        matrix.char_to_index = char_to_index
        matrix.values = values
        return matrix
      
    @classmethod
    def create_from_other( Class, other, values=None ):
        """
        Create a new Matrix with attributes taken from `other` but with the 
        values taken from `values` if provided
        """
        m = Class()
        m.alphabet = other.alphabet
        m.sorted_alphabet = other.sorted_alphabet
        m.char_to_index = other.char_to_index
        if values is not None:
            m.values = values
        else:
            m.values = other.values
        return m
            
    @property
    def width( self ):
        """
        Return the width (size along the sequence axis) of this matrix.
        """
        return self.values.shape[0]

    def reverse_complement( self ):
        """
        Create the reverse complement of this matrix. The result probably
        only makese sense if the alphabet is that of DNA ('A','C','G','T').
        """
        rval = copy( self )
        # Conveniently enough, reversing rows and columns is exactly what we
        # want, since this results in A swapping with T and C swapping with G.
        rval.values = self.values[::-1,::-1].copy()
        return rval

class FrequencyMatrix( BaseMatrix ):
    """
    A position specific count/frequency matrix.
    """
        
    DEFAULT_CORRECTION = 0.0000000001
    """
    Default value to use for correcting when dealing with counts of zero,
    chosen to produce scoring matrices that are the same as produced by CREAD.
    """
        
    def to_logodds_scoring_matrix( self, background=None, correction=DEFAULT_CORRECTION ):
        """
        Create a standard logodds scoring matrix.
        """
        alphabet_size = len( self.alphabet )
        if background is None:
            background = ones( alphabet_size, float32 ) / alphabet_size
        # Row totals as a one column array
        totals = sum( self.values, 1 )[:,newaxis]
        values = log2( maximum( self.values, correction ) ) \
               - log2( totals ) \
               - log2( maximum( background, correction ) )
        return ScoringMatrix.create_from_other( self, values.astype( float32 ) )

    def to_stormo_scoring_matrix( self, background=None ):
        """
        Create a scoring matrix from this count matrix using the method from:

        Hertz, G.Z. and G.D. Stormo (1999). Identifying DNA and protein patterns with statistically 
        significant alignments of multiple sequences. Bioinformatics 15(7): 563-577.
        """
        alphabet_size = len( self.alphabet )
        if background is None:
            background = ones( alphabet_size, float32 ) / alphabet_size
        # Row totals as a one column array
        totals = sum( self.values, 1 )[:,newaxis]
        values = log2( self.values + background ) \
               - log2( totals + 1 ) - log2( background )
        return ScoringMatrix.create_from_other( self, values.astype( float32 ) )
        
class ScoringMatrix( BaseMatrix ):
    """
    A position specific matrix containing values that are suitable for
    scoring a sequence.
    """
        
    def score_string( self, string ):
        """
        Score each valid position in `string` using this scoring matrix. 
        Positions which were not scored are set to nan.
        """
        rval = zeros( len( string ), float32 )
        rval[:] = nan
        _pwm.score_string( self.values, self.char_to_index, string, rval )
        return rval
        
    def score_string_with_gaps( self, string ):
        """
        Score each valid position in `string` using this scoring matrix. 
        Positions which were not scored are set to nan. Gap characters are
        ignored (matrices score across them).
        """
        rval = zeros( len( string ), float32 )
        rval[:] = nan
        _pwm.score_string_with_gaps( self.values, self.char_to_index, string, rval )
        return rval