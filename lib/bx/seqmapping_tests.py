"""
Tests for `bx.seqmapping`.
"""

import sys

import unittest
import bx.seqmapping

from numpy import array, allclose
from StringIO import StringIO

class CharMappingTests( unittest.TestCase ):
    def test_DNA( self ):
        assert( allclose( bx.seqmapping.DNA.translate( "ACGTacgt-?X" ),
                          [ 0, 1, 2, 3, 0, 1, 2, 3, 4, -1, -1 ] ) )
    def test_DNA_list( self ):
        assert( allclose( bx.seqmapping.DNA.translate_list( [ "ACGTA", "TGCAX" ] ),
                          [ 0 + 3*6, 1 + 2*6, 2 + 1*6, 3 + 0*6, -1 ] ) )
    def test_other( self ):
        m = bx.seqmapping.CharToIntArrayMapping()
        m.set_mapping( "A", 0 )
        m.set_mapping( "B", 7 )
        assert( allclose( m.translate( "ABCCBA" ), [ 0, 7, -1, -1, 7, 0 ] ) )
        
class IntMappingTests( unittest.TestCase ):
    def test_simple( self ):
        m = bx.seqmapping.IntToIntMapping( 4 )
        m.set_mapping( 0, 0 )
        m.set_mapping( 2, 0 )
        m.set_mapping( 1, 1 )
        m.set_mapping( 3, 1 )
        assert( allclose( m.translate( array( [ 0, 1, 2, 3, 4 ], 'i' ) ), array( [ 0, 1, 0, 1, -1 ] ) ) )
   
eight_species_mapping = """TTTTTTTT 0
CCCCCCCC 4
AAAAAAAA 1
GGGGGGGG 5
AAAAA*AA 2
TTTTT*TT 0
GGGGG*GG 5
CCCCC*CC 4
GGGG*GGG 5
TTTT*TTT 2
GGGAAGGG 5
AAAA*AAA 2
TTTTT*T* 2
CCCCC*C* 4
CCCTTCCC 4
CCCC*CCC 4
TTTT**TT 2
AAAA**AA 2
AAAAA*A* 2
GGGGG*G* 5
AAAAAGAA 2
TTTTTCTT 2
GGGAA*GG 5
TTTT**T* 2
TTTCCTTT 0
AAAAAAA* 1
CCCTT*CC 3
TTTTTTT* 2
CC*CC*CC 3
AAAGGAAA 2
------G- 2
"""
       
rows = [ "AAATTGT-----ATGTCCATCCTTTAAAGGTCATTCCTTTAATGGTCTTTTCTGGACACCACTAGGGGTCAGAAGTAGTTCATCAAAC-----------------TTTCTTCCCTCCC-TACTTCAGTG",
         "AAATTGT-----ATGTCCATCCTTTAAAGGTCATTCCTTTAATGGTCTTTTCTGGACACCACTAGGGGTCAGAAGTAGTTCATCAAAC-----------------TTTCTTCCCTCCC-TACTTCAGTG",
         "AAATTTT-----ATGTCTATCCTTTAAAGGTCATTCCTCTAATAGTCTTTTCTGGACACCACTAGGGGTCAGAAGTAGTTCATTAAAC-----------------TTTCTTCCCTCCC-TACCTCAGTG",
         "AAACTGT-----ATCACCACCTTTTTAAGGTCATTTCTCTAATGATCCTGTT-GCATACCAGTAGGGGGCAGAAGTGTTCCGCTGATTTCCGCCCTCCTCCCCACCCCCCCACCCCCC-TTATTCAAAG",
         "*********************************************************************************************************************************",
         "-TATTAT-----ATGGCCATGTTCAAAAGGTTGTTTCTCTAATGATTCCTTC-TGATACCAGTAGGGGTCAGAAGTGGTCCATTGATT---------------------CTTTTCCTC-TGATTC-AAG",
         "AAATTGA--AAGATCTCACTCTTTGCCAGGTAGTCCATCTAAGGGTCACATATGGATACCAGCAGGGCCT-GAAGAAGCCCATTGAAT------------------------TTTCCC-ATCTTCAAGG",
         "AAATTCATGATAGTGTCACTCTTAAATAGATGATTC--------TTCACAT---GATGCCAGCAGGGGGC-AGAGCAGGCTGTGAAAT------------------------TTTCCCTTTCTTCAAAG" ]

class AlignmentMappingTests( unittest.TestCase ):
    def test_largescale( self ):
       f = StringIO( eight_species_mapping )
       n, m = bx.seqmapping.alignment_mapping_from_file( f )
       t = bx.seqmapping.DNA.translate_list( rows )
       i = m.translate( t )
        
        
test_classes = [ AlignmentMappingTests, CharMappingTests, IntMappingTests ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
