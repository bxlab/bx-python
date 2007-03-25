"""
Numeric arrays stored as individually compressed blocks on disk, allowing
pseudo-random acccess. 

`BinnedArray` is used to build such an array in memory and save it to disk.
`BinnedArrayWriter` can instead be used when creating the array sequentially
(does not require keeping all data in memory). `FileBinnedArray` provides
read only access to an on disk binned array. 
"""

from __future__ import division

import math

from numpy import *
from struct import *
from bx_extras.lrucache import LRUCache

import sys
platform_is_little_endian = ( sys.byteorder == 'little' )

MAGIC=0x4AB04612

# Version incremented from version 0 to version 1 by Ian Schenck, June
# 23, 2006.  Version 1 supports different typecodes, and in doing so
# breaks the original header format.  The new FileBinnedArray is
# backwards compatible with version 0.

# Version 1 -> 2 by James Taylor, allow specifying different compression 
# types.

VERSION=2

# Compression types

comp_types = dict()

comp_types['none'] = ( lambda x: x, lambda x: x )

try:
    import zlib
    comp_types['zlib'] = ( zlib.compress, zlib.decompress )
except:
    pass
    
try:
    import lzo
    comp_types['lzo'] = ( lzo.compress, lzo.decompress )
except:
    pass


MAX=512*1024*1024 

class BinnedArray( object ):
    def __init__( self, bin_size=512*1024, default=NaN, max_size=MAX, typecode="f" ):
        self.max_size = max_size
        self.bin_size = bin_size
        self.nbins = int( math.ceil( ( max_size / self.bin_size ) ) )
        self.bins = [ None ] * self.nbins
        self.default = default
        self.typecode = typecode
    def get_bin_offset( self, index ):
        return index // self.bin_size, index % self.bin_size
    def init_bin( self, index ):
        # self.bins[index] = zeros( self.bin_size ) * self.default
        self.bins[index] = zeros( self.bin_size, self.typecode )
        self.bins[index][:] = self.default
    def get( self, key ):
        bin, offset = self.get_bin_offset( key )
        if self.bins[bin] is None:
            return self.default
        else:
            return self.bins[bin][offset]
    def set( self, key, value ):
        bin, offset = self.get_bin_offset( key )
        if self.bins[bin] is None: 
            self.init_bin( bin )
        self.bins[bin][offset] = value
    def get_range( self, start, end ):
        size = end - start
        assert size >= 0
        rval = []
        while size > 0:
            bin, offset = self.get_bin_offset( start )
            delta = self.bin_size - offset
            if self.bins[bin] is None:
                if delta < size:
                    rval.append( resize( array(self.default, self.typecode), (delta,) ) )
                    size -= delta
                    start += delta
                else:
                    rval.append( resize( array(self.default, "f"), (size,) ) )
                    size = 0
            else:
                if delta < size:
                    rval.append( self.bins[bin][offset:offset+delta] )
                    size -= delta
                    start += delta
                else:
                    rval.append( self.bins[bin][offset:offset+size] )
                    size = 0
        return concatenate( rval )
    def __getitem__( self, key ):
        if isinstance( key, slice ):
            start, stop, stride = key.indices( self.max_size )
            assert stride == 1, "Slices with strides are not supported"
            return self.get_range( start, stop )
        else:
            return self.get( key )
    def __setitem__( self, key, value ):
        return self.set( key, value )
    def to_file( self, f, comp_type='zlib' ):
        # Get compress method
        compress, _ = comp_types[comp_type]
        # Write header
        write_packed( f, ">5I", MAGIC, VERSION, self.max_size, self.bin_size, self.nbins )
        # save type code
        f.write( pack('c',self.typecode ) )
        # save compression type
        f.write( comp_type[0:4].ljust( 4 ) )
        # write default value
        a = array( self.default, self.typecode ) 
        # Struct module can't deal with NaN and endian conversion, we'll hack 
        # around that by byteswapping the array
        if platform_is_little_endian: 
            a = a.byteswap()
        f.write( a.tostring() )
        # Save current position (start of bin offsets)
        index_start_pos = f.tell()
        # Skip forward to save space for index
        f.seek( calcsize( ">2I" ) * self.nbins, 1 )
        bin_pos_and_size = []
        # Write each bin
        for bin in self.bins:
            if bin is None: 
                bin_pos_and_size.append( ( 0, 0 ) )
            else:
                assert bin.dtype.char == self.typecode
                if platform_is_little_endian:
                    s = bin.byteswap().tostring()
                else:
                    s = bin.tostring()
                compressed = compress( s )
                bin_pos_and_size.append( ( f.tell(), len( compressed ) ) )
                f.write( compressed )
        # Go back and fill in table
        f.seek( index_start_pos )
        for pos, size in bin_pos_and_size:
            write_packed( f, ">2I", pos, size )
            
class FileBinnedArray( object ):
    def __init__( self, f, cache=32):
        # If cache=None, then everything is allowed to stay in memory,
        # this is the default behavior.
        self.f = f
        M, V, max_size, bin_size, nbins = read_packed( f, ">5I" )
        assert M == MAGIC
        # assert version less than max supported
        assert V <= VERSION, "File is version %d but I don't know about anything beyond %d" % ( V, VERSION )
        self.max_size = max_size
        self.bin_size = bin_size
        self.nbins = nbins        
        self.bins = LRUCache(size=cache)
        # Read typecode
        if V >= 1:
            self.typecode = unpack( 'c', f.read(1) )[0]
        else:
            self.typecode = 'f'
        # Read compression type
        if V >= 2:
            self.comp_type = f.read( 4 ).strip()
        else:
            self.comp_type = 'zlib'
        self.decompress = comp_types[self.comp_type][1]
        # Read default value
        s = f.read( calcsize( self.typecode ) )
        a = fromstring( s, self.typecode )
        if platform_is_little_endian: 
            a = a.byteswap()
        self.default = a[0]
        # Read bin sizes and offsets
        self.bin_pos = []
        self.bin_sizes = []
        for i in range( nbins ):
            pos, size = read_packed( f, ">2I" )
            self.bin_pos.append( pos )
            self.bin_sizes.append( size )
    def get_bin_offset( self, index ):
        return int( index // self.bin_size ), int( index % self.bin_size )
    def load_bin( self, index ):
        assert self.bin_pos[index] != 0
        self.f.seek( self.bin_pos[index] )
        raw = self.f.read( self.bin_sizes[index] )
        a = fromstring( self.decompress( raw ), self.typecode )
        if platform_is_little_endian:
            a = a.byteswap()
        assert len( a ) == self.bin_size
        self.bins[index] = a
    def get( self, key ):
        bin, offset = self.get_bin_offset( key )
        if bin in self.bins:
            return self.bins[bin][offset]
        elif self.bin_pos[bin]:
            self.load_bin( bin )
            return self.bins[bin][offset]
        else:
            return self.default
    def get_range( self, start, end ):
        size = end - start
        assert size >= 0
        rval = []
        while size > 0:
            bin, offset = self.get_bin_offset( start )
            delta = self.bin_size - offset
            if not bin in self.bins and self.bin_pos[bin] != 0:
                self.load_bin( bin )
            if self.bins[bin] is None:
                if delta < size:
                    rval.append( resize( array(self.default, self.typecode), (delta,) ) )
                    size -= delta
                    start += delta
                else:
                    rval.append( resize( array(self.default, self.typecode), (size,) ) )
                    size = 0
            else:
                if delta < size:
                    rval.append( self.bins[bin][offset:offset+delta] )
                    size -= delta
                    start += delta
                else:
                    rval.append( self.bins[bin][offset:offset+size] )
                    size = 0
        return concatenate( rval )
    def __getitem__( self, key ):
        if isinstance( key, slice ):
            start, stop, stride = key.indices( self.max_size )
            assert stride == 1, "Slices with strides are not supported"
            return self.get_range( start, stop )
        else:
            return self.get( key )
        
class BinnedArrayWriter( object ):
    def __init__( self, f, bin_size=512*1024, default=NaN, max_size=MAX, typecode="f", comp_type='zlib' ):
        # All parameters in the constructor are immutable after creation
        self.f = f
        self.max_size = max_size
        self.bin_size = bin_size
        self.nbins = int( math.ceil( ( max_size / self.bin_size ) ) )
        self.default = default
        self.typecode = typecode
        self.bin = 0
        self.bin_pos = 0
        self.bin_index = []
        self.buffer = resize( array(self.default, self.typecode), (self.bin_size,) )
        self.buffer_contains_values = False
        self.comp_type = comp_type
        self.compress = comp_types[comp_type][0]
        self.write_header()
        # Start the first bin
        ## self.bin_index = [ (self.data_offset, 0) ]
        # Put the fp at the start of the data (we go back and fill in the index at the end)
        self.f.seek( self.data_offset )

    def write_header( self ):
        self.f.seek(0)
        # Write header
        write_packed( self.f, ">5I", MAGIC, VERSION, self.max_size, self.bin_size, self.nbins )
        # save type code
        self.f.write( pack('c',self.typecode ) )
        # write default value
        a = array( self.default, self.typecode ) 
        # write comp type
        self.f.write( self.comp_type[0:4].ljust(4) )
        # write default
        # Struct module can't deal with NaN and endian conversion, we'll hack 
        # around that by byteswapping the array
        if platform_is_little_endian: 
            a = a.byteswap()
        self.f.write( a.tostring() )
        # Save current position (start of bin offsets)
        self.index_pos = self.f.tell()
        self.data_offset = self.index_pos + (self.nbins * calcsize( ">2I" ))
        
    def write_index( self ):
        self.f.seek(self.index_pos)
        for pos, size in self.bin_index:
            write_packed( self.f, ">2I", pos, size )

    def skip( self ):
        self.bin_pos += 1
        if self.bin_pos == self.bin_size:
            self.flush()
            self.bin_pos = 0
            self.bin += 1
            assert self.bin <= self.nbins
            self.buffer = resize( array(self.default, self.typecode), (self.bin_size,) )
            self.buffer_contains_values = False
            ## self.bin_index.append( (self.f.tell(), 0) )

    def write( self, data ):
        self.buffer[self.bin_pos] = data
        self.buffer_contains_values = True
        self.bin_pos += 1
        if self.bin_pos == self.bin_size:
            self.flush()
            self.bin_pos = 0
            self.bin += 1
            assert self.bin <= self.nbins
            self.buffer = resize( array(self.default, self.typecode), (self.bin_size,) )
            self.buffer_contains_values = False
            ## self.bin_index.append( (self.f.tell(), 0) )

    def flush( self ):
        # Flush buffer to file
        if self.buffer_contains_values:
            ## pos, size = self.bin_index[self.bin]
            ## self.f.seek( pos )
            pos = self.f.tell()
            if platform_is_little_endian:
                s = self.buffer.byteswap().tostring()
            else:
                s = self.buffer.tostring()
            compressed = self.compress( s )
            size = len( compressed )
            assert len( self.bin_index ) == self.bin
            self.bin_index.append( ( pos, size ) )
            self.f.write( compressed )
        else:
            assert len( self.bin_index ) == self.bin
            self.bin_index.append( ( 0, 0 ) )

    def finish( self ):
        self.flush()
        self.nbins = self.bin + 1
        self.write_header()
        self.write_index()

def write_packed( f, pattern, *vals ):
    f.write( pack( pattern, *vals ) )
    
def read_packed( f, pattern ):
    rval = unpack( pattern, f.read( calcsize( pattern ) ) )
    if len( rval ) == 1: return rval[0]
    return rval