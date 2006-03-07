"""
Functions for reading and writing 2d matrices in the format used
by SVDPACK / SVDLIBC and also easily parsed by R
"""

from numarray import *

import struct

header_pattern = "!ll" 
header_pattern_size = struct.calcsize( header_pattern )

def write_dense_binary( a, f ):
    """Write the array `a` to the file `f` in dense binary format"""
    assert len( a.shape ) == 2, "Only 2 dimensional arrays allowed"
    nrow, ncol = a.shape
    # Write shape first as two longs
    f.write( struct.pack( header_pattern, nrow, ncol) )    
    # Write each entry as a float. array.tofile does not allow control of byte order
    row_pat = "!%df" % ncol
    row_pat_size = struct.calcsize( row_pat )
    for i in range( nrow ):
        f.write( struct.pack( row_pat, *a[i] ) )
   
def read_dense_binary( f ):
    """Read a 2d Float32 array from file `f` and return it"""
    nrow, ncol = struct.unpack( header_pattern, f.read( header_pattern_size ) )
    a = zeros( ( nrow, ncol ), type=Float32 )
    nbytes = nrow * ncol * 4 
    buf  = f.read( nbytes )
    if len( buf ) != nbytes: raise "Didn't read expected number of bytes" 
    return fromstring( buf, Float32, ( nrow, ncol ) )    
    #row_pat = "!%df" % ncol
    #row_pat_size = struct.calcsize( row_pat )
    #for i in range( nrow ):
    #    a[i] = array( struct.unpack( row_pat, f.read( row_pat_size ) ) )
    #return a

class DBRowReader( object ):
    """Read a 2d Float32 array from file `f` and return it row by row"""
    def __init__( self, f ):
        self.f = f
        self.nrow, self.ncol = struct.unpack( header_pattern, f.read( header_pattern_size ) )
        self.row_pat = "!%df" % self.ncol
        self.row_pat_size = struct.calcsize( self.row_pat )
        self.rows_read = 0
    def close( self ):
        self.f.close()
    def next( self ):
        if self.rows_read >= self.nrow: raise StopIteration
        self.rows_read += 1
        return array( struct.unpack( self.row_pat, self.f.read( self.row_pat_size ) ) )
    def __iter__( self ): return self

class DBRowWriter( object ):
    def __init__( self, nrow, ncol, f ):
        self.nrow = nrow
        self.ncol = ncol
        self.f = f
        self.row_pat = "!%df" % self.ncol
        self.row_pat_size = struct.calcsize( self.row_pat )
        self.rows_written = 0
        # Write shape first as two longs
        self.f.write( struct.pack( header_pattern, self.nrow, self.ncol) )    
    def write( self, row ):
        assert len( row ) == self.ncol, "Too many columns in row"
        assert self.rows_written < self.nrow, "Too many rows written"
        # Write each entry as a float. array.tofile does not allow control of byte order
        self.f.write( struct.pack( self.row_pat, *row ) )
        self.rows_written += 1
    def close( self ):
        assert self.rows_written == self.nrow
        self.f.close()

def test():
    import numarray.random_array, StringIO
    f = StringIO.StringIO()
    a = numarray.random_array.random( ( 500, 300 ) )
    write_dense_binary( a, f )
    f.seek( 0 )
    b = read_dense_binary( f )
    assert allclose( a, b )
    


if __name__ == "__main__": test()
