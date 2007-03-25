"""
Classes for char-to-int mapping and int-to-int mapping.

:Author: James Taylor (james@bx.psu.edu)

The char-to-int mapping can be used to translate a list of strings
over some alphabet to a single int array (example for encoding a multiple
sequence alignment).

The int-to-int mapping is particularly useful for creating partitions,
and provides methods to merge/split symbols in the output mapping.

The two forms of mapping can be combined, for example to encode a 
multiple sequence alignment in a reduced alphabet defined by a partition
of alignment columns. Many of the helper methods provided are for 
solving such alignment oriented problems. 

This code was originally written for the `ESPERR`_ project which includes
software for searcing for alignment encodings that work well for specific
classification problems using various Markov chain classifiers over the
reduced encodings.

Most of the core implementation is in the pyrex/C extension 
"_seqmapping.pyx" for performance reasons (specifically to avoid the
excessive bounds checking that would make a sequence/array lookup heavy
problem like this slow in pure python).

.. _ESPERR: http://www.bx.psu.edu/projects/esperr/
"""

from _seqmapping import *

# Char->Int mapping for DNA characters with missing data
                
DNA = CharToIntArrayMapping()
DNA.set_mapping( "a", 0 )
DNA.set_mapping( "A", 0 )
DNA.set_mapping( "c", 1 )
DNA.set_mapping( "C", 1 )
DNA.set_mapping( "g", 2 )
DNA.set_mapping( "G", 2 )
DNA.set_mapping( "t", 3 )
DNA.set_mapping( "T", 3 )
DNA.set_mapping( "-", 4 )
DNA.set_mapping( "*", 5 )

# Creating mappings

def alignment_mapping_from_file( f, char_mapping=DNA ):
    """
    Create a mapping from a file of alignment columns.
    """    
    columns, symbols = [], []
    for line in f:
        column, symbol = line.split()
        columns.append( column )
        symbols.append( int( symbol ) )
                
    align_count = len( columns[0] )
        
    mapping = IntToIntMapping( char_mapping.get_out_size() ** align_count )
        
    for column, symbol in zip( columns, symbols ):
        index = char_mapping.translate_list( list( column ) )[0]
        mapping.set_mapping( index, symbol )

    return align_count, mapping

def second_mapping_from_file( f, first_mapping, char_mapping=DNA ):
        
    columns, symbols = [], []
    for line in f:
        column, symbol = line.split()
        columns.append( column )
        symbols.append( int( symbol ) )
                
    align_count = len( columns[0] )
        
    mapping = IntToIntMapping( first_mapping.get_out_size() )
        
    for column, symbol in zip( columns, symbols ):
        index = char_mapping.translate_list( list( column ) )[0]
        if first_mapping[index] >= 0:
            mapping.set_mapping( first_mapping[index], symbol )

    return mapping


def identity_mapping( size ):
    mapping = IntToIntMapping( size )
    for i in range( size ):
        mapping.set_mapping( i, i )
    return mapping    
    
