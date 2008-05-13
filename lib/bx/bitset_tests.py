"""
Tests for `bx.bitset`.
"""

import bx.bitset
import unittest

class AbstractTests( object ):

    def assert_bits( self, bits, list ):
        assert bits.size == len( list ), "Bitset size and verification list size do not match"
        for i in range( bits.size ):
            self.assertEquals( bits[i], list[i] )

    def test_overflow_create( self ):
        self.assertRaises( ValueError, self.new_bits, 4000000000 )
        
    def test_overflow_access( self ):
        bits = self.new_bits( 100 )
        self.assertRaises( IndexError, bits.set, -5 )
        self.assertRaises( IndexError, bits.set, 110 )

    def test_access( self ):
        # Create and assert empty
        bits = self.new_bits( 100 )
        l = [ 0 ] * 100
        self.assert_bits( bits, l )
        # Set some positions
        for pos in ( 11, 14, 70, 16 ):
            bits.set( pos )
            l[ pos ] = 1
        # Clear some positions
        for pos in ( 14, 80, 16 ):
            bits.clear( pos )
            l[ pos ] = 0
        self.assert_bits( bits, l )

    def test_range_access( self ):
        # Create and assert empty
        bits = self.new_bits( 100 )
        l = [ 0 ] * 100
        self.assert_bits( bits, l )
        # Set some positions
        for b, e in ( ( 11, 14 ), (20,75), (90,99) ):
            bits.set_range( b, e-b)
            for pos in range( b, e ): l[ pos ] = 1
        self.assert_bits( bits, l )

    def test_count( self ):
        # Create and assert empty
        bits = self.new_bits( 100 )
        # Set some positions
        for b, e in ( ( 11, 14 ), (20,75), (90,100) ):
            bits.set_range( b, e-b)
        self.assertEquals( bits.count_range( 0, 0 ), 0 )
        self.assertEquals( bits.count_range( 0, 20 ), 3 )
        self.assertEquals( bits.count_range( 25, 25 ), 25 )
        self.assertEquals( bits.count_range( 80, 20 ), 10 )
        self.assertEquals( bits.count_range( 0, 100 ), 68 )

    def test_find( self ):
        # Create and assert empty
        bits = self.new_bits( 100 )
        # Set some positions
        for b, e in ( ( 11, 14 ), (20,75), (90,100) ):
            bits.set_range( b, e-b)
        # Next set
        self.assertEquals( bits.next_set( 0 ), 11 )
        self.assertEquals( bits.next_set( 13 ), 13 )
        self.assertEquals( bits.next_set( 15 ), 20 )
        # Next clear
        self.assertEquals( bits.next_clear( 0 ), 0 )
        self.assertEquals( bits.next_clear( 11 ), 14 )
        self.assertEquals( bits.next_clear( 20 ), 75 )
        self.assertEquals( bits.next_clear( 92 ), 100 )

    def test_and( self ):
        bits1 = self.new_bits( 100 )
        bits2 = self.new_bits( 100 )
        bits1.set_range( 20, 40 )
        bits2.set_range( 50, 25 )
        bits1.iand( bits2 )
        l = [0]*100
        for i in range( 50, 60 ): l[i] = 1
        self.assert_bits( bits1, l )
    
    def test_or( self ):
        bits1 = self.new_bits( 100 )
        bits2 = self.new_bits( 100 )
        bits1.set_range( 20, 40 )
        bits2.set_range( 50, 25 )
        bits1.ior( bits2 )
        l = [0]*100
        for i in range( 20, 75 ): l[i] = 1
        self.assert_bits( bits1, l )
        
    def test_not( self ):
        bits = self.new_bits( 100 )
        bits.set_range( 20, 40 )
        bits.invert()
        l = [1]*100
        for i in range( 20, 60 ): l[i] = 0
        self.assert_bits( bits, l )
        
class BitSetTests( AbstractTests, unittest.TestCase ):
    def new_bits( self, size ):
        return bx.bitset.BitSet( size ) 

class BinnedBitSetTests( AbstractTests, unittest.TestCase ):
    def new_bits( self, size ):
        granularity = size % 11 
        return bx.bitset.BinnedBitSet( size, granularity ) 