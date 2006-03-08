"""
Unit tests for seqmapping.py
"""

import unittest
import bx.seqmapping

from Numeric import array

class CharMappingTests( unittest.TestCase ):
    def test_DNA( self ):
        self.assertEqual( bx.seqmapping.DNA.translate( "ACGTacgt-?X" ),
                          [ 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, -1 ] )
    def test_DNA_list( self ):
        self.assertEqual( bx.seqmapping.DNA.translate_list( [ "ACGTA", "TGCAX" ] ),
                          [ 0 + 3*6, 1 + 2*6, 2 + 1*6, 3 + 0*6, -1 ] )
    def test_other( self ):
        m = bx.seqmapping.CharToIntArrayMapping()
        m.set_mapping( "A", 0 )
        m.set_mapping( "B", 7 )
        self.assertEqual( m.translate( "ABCBCA" ), [ 0, 7, -1, -1, 7, 0 ] )
        
class IntMappingTests( unittest.TestCase ):
    def test_simple( self ):
        m = bx.seqmapping.IntToIntMapping( 4 )
        m.set_mapping( 0, 0 )
        m.set_mapping( 2, 0 )
        m.set_mapping( 1, 1 )
        m.set_mapping( 3, 1 )
        self.assertEqual( m.translate( array( [ 0, 1, 2, 3, 4 ], 'i' ) ), array( [ 0, 1, 0, 1, -1 ] ) ) 
        
test_classes = [ CharMappingTests, IntMappingTests ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
