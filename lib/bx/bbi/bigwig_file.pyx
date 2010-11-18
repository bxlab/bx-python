"""
BigWig file.
"""

from bbi_file cimport BBIFile

DEF big_wig_sig = 0x888FFC26

cdef class BigWigFile( BBIFile ): 
    """
    A "big binary indexed" file whose raw data is in wiggle format.
    """
    def __init__( self, file=None ):
        BBIFile.__init__( self, file, big_wig_sig, "bigwig" )
        