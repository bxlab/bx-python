"""
Access to files containing sequence data in 'twobit' format.
"""

import sys
import _twobit

from struct import *
from UserDict import DictMixin

TWOBIT_MAGIC_NUMBER = 0x1A412743
TWOBIT_MAGIC_NUMBER_SWAP = 0x4327411A
TWOBIT_MAGIC_SIZE = 4

TWOBIT_VERSION = 0

class TwoBitSequence( object ):
    def __init__( self, tbf, header_offset=None ):
        self.tbf = tbf
        self.header_offset = header_offset
        self.sequence_offset = None
        self.size = None
        self.n_blocks = None
        self.masked_blocks = None
        self.loaded = False
        
    def __getitem__( self, slice ):
        start, stop, stride = slice.indices( self.size )
        assert stride == 1, "Striding in slices not supported"
        if stop - start < 1:
            return ""
        return _twobit.read( self.tbf.file, self, start, stop, self.tbf.do_mask  )
        
    def __len__( self ):
        return self.size
        
    def get( self, start, end ):
        # Trim start / stop
        if start < 0:
            start = 0
        if end > self.size:
            end = self.size
        out_size = end - start
        if out_size < 1:
            raise Exception( "end before start (%s,%s)" % ( start,end ) )
        # Find position of packed portion
        dna = _twobit.read( self.tbf.file, self, start, end, self.tbf.do_mask )
        # Return
        return dna
        
class TwoBitFile( DictMixin ):
    def __init__( self, file, do_mask=True ):
        self.do_mask = do_mask
        # Read magic and determine byte order
        self.byte_order = ">"
        magic = unpack( ">L", file.read( TWOBIT_MAGIC_SIZE ) )[0]
        if magic != TWOBIT_MAGIC_NUMBER:
            if magic == TWOBIT_MAGIC_NUMBER_SWAP: 
                self.byte_order = "<"
            else: 
                raise Exception( "Not a NIB file" )
        self.magic = magic
        self.file = file
        # Read version
        self.version = self.read( "L" )
        if self.version != TWOBIT_VERSION:
            raise Exception( "File is version '%d' but I only know about '%d'" % ( self.version, TWOBIT_VERSION ) )
        # Number of sequences in file
        self.seq_count = self.read( "L" )
        # Header contains some reserved space
        self.reserved = self.read( "L" )
        # Read index of sequence names to offsets
        index = dict()
        for i in range( self.seq_count ):
            name = self.read_p_string()
            offset = self.read( "L" )
            index[name] = TwoBitSequence( self, offset )
        self.index = index

    def __getitem__( self, name ):
        seq = self.index[name]
        if not seq.loaded:
            self.load_sequence( name )
        return seq
        
    def keys( self ):
        return self.index.keys()
        
    def load_sequence( self, name ):
        seq = self.index[name]
        # Seek to start of sequence block
        self.file.seek( seq.header_offset )
        # Size of sequence
        seq.size = self.read( "L" )
        # Read N and masked block regions
        seq.n_block_starts, seq.n_block_sizes = self.read_block_coords()
        seq.masked_block_starts, seq.masked_block_sizes = self.read_block_coords()
        # Reserved
        self.read( "L" )
        # Save start of actualt sequence
        seq.sequence_offset = self.file.tell()
        # Mark as loaded
        seq.loaded = True
        
    def read_block_coords( self ):
        blocks = []
        block_count = self.read( "L" )
        if block_count == 0:
            return [], []
        starts = self.read( str( block_count ) + "L", untuple=False )
        sizes = self.read( str( block_count ) + "L", untuple=False  )
        return list( starts ), list( sizes )
        
    def read( self, pattern, untuple=True ):
        rval = unpack( self.byte_order + pattern, 
                       self.file.read( calcsize( self.byte_order + pattern ) ) )
        if untuple and len( rval ) == 1: 
            return rval[0]
        return rval
        
    def read_p_string( self ):
        """
        Read a length-prefixed string 
        """
        length = self.read( "B" )
        return self.file.read( length )
