from __future__ import division

import sys, struct, string, math

NIB_MAGIC_NUMBER = 0x6BE93D3A
NIB_MAGIC_NUMBER_SWAP = 0x3A3DE96B
NIB_MAGIC_SIZE = 4
NIB_LENGTH_SIZE = 4
NIB_I2C_TABLE = "TCAGNXXXtcagnxxx"

READ_CHUNK_SIZE=1024*1024

class NibFile( object ):
    def __init__( self, file ):
        self.byte_order = ">" 
        magic = struct.unpack( ">L", file.read( NIB_MAGIC_SIZE ) )[0]
        if magic != NIB_MAGIC_NUMBER: 
            if magic == NIB_MAGIC_NUMBER_SWAP: self.byte_order = "<"  
            else: raise "Not a NIB file"
        self.file = file
        self.magic = magic
        self.length = struct.unpack( "%sL" % self.byte_order, file.read( NIB_LENGTH_SIZE ) )[0]
    def get( self, start, length ):
        assert start >= 0
        assert start + length - 1 < self.length
        # Read block of bytes containing sequence
        block_start = int( math.floor( start / 2 ) )
        block_end = int( math.floor( ( start + length - 1 ) / 2 ) )
        block_len = block_end + 1 - block_start
        self.file.seek( NIB_MAGIC_SIZE + NIB_LENGTH_SIZE + block_start )
        result = []
        raw = self.file.read( block_len )
        data = struct.unpack( "%s%dB" % ( self.byte_order, block_len ), raw )
        # Translate to character representation
        for value in data:
            result.append( NIB_I2C_TABLE[ ( value >> 4 ) & 0xF ] )
            result.append( NIB_I2C_TABLE[ ( value >> 0 ) & 0xF ] )
        # Trim if start / end are odd 
        if start & 1: del result[ 0 ]
        if ( start + length ) & 1: del result[ -1 ]
        # Return as string
        return string.join( result, '' )
