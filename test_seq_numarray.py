from numarray import *
import unittest
from bx import seq_numarray

class TestCase( unittest.TestCase ):
    def testBasic( self ):
        assert array_equal( array( [ 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 ] ),
                            seq_numarray.DNA.translate( "-ACGT-ACGT" ) )
    def testCol( self ):
        assert seq_numarray.DNA.translate_alignment_column( "---" ) == 0
        assert seq_numarray.DNA.translate_alignment_column( "TTT" ) == 124
        
    def testReverseCol( self ):
        assert seq_numarray.DNA.reverse_alignment_column( 3, 0 ) == "---"
        assert seq_numarray.DNA.reverse_alignment_column( 3, 124 ) == "TTT"

if __name__ == "__main__":
    unittest.main()
