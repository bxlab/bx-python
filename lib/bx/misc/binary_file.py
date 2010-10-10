"""
Wrappers for doing binary IO on file-like objects
"""

import numpy
import struct
import sys

## Standard size:
## short is 8 bits
## int and long are 32 bits
## long long is 64 bits

class BadMagicNumber( IOError ):
    pass

class BinaryFileReader( object ):
    """
    Wrapper for doing binary reads on any file like object.
    
    Currently this is not heavily optimized (it uses the `struct` module to
    unpack)
    """    
    def __init__( self, file, magic = None, is_little_endian = False ):
        self.is_little_endian = is_little_endian
        self.file = file
        if magic is not None:
            # Attempt to read magic number and chuck endianess
            bytes = file.read( 4 )
            if struct.unpack( ">I", bytes )[0] == magic:
                pass
            elif struct.unpack( "<I", bytes )[0] == magic:
                self.is_little_endian = True
            else:
                raise BadMagicNumber( "File does not have expected magic number: %x != %x or %x" \
                        % ( magic, struct.unpack( ">I", bytes )[0], struct.unpack( "<I", bytes )[0] ) )
        # Set endian code
        if self.is_little_endian:
            self.endian_code = "<"
            self.byteswap_needed = ( sys.byteorder != "little" )
        else:
            self.endian_code = ">"
            self.byteswap_needed = ( sys.byteorder != "big" )
        
    def unpack( self, format, buffer, byte_count=None ):
        pattern = "%s%s" % ( self.endian_code, format )
        if byte_count is None:
            byte_count = struct.calcsize( pattern )
        return struct.unpack( pattern, buffer )
        
    def read_and_unpack( self, format, byte_count=None ):
        """
        Read enough bytes to unpack according to `format` and return the
        tuple of unpacked values.
        """
        pattern = "%s%s" % ( self.endian_code, format )
        if byte_count is None:
            byte_count = struct.calcsize( pattern )
        return struct.unpack( pattern, self.file.read( byte_count ) )
    
    def read_c_string( self ):
        """
        Read a zero terminated (C style) string
        """
        rval = []
        while 1:
            ch = self.file.read(1)
            assert len( ch ) == 1, "Unexpected end of file"
            if ch == '\0':
                break
            rval.append( ch )
        return ''.join( rval )
        
    def read_raw_array( self, dtype, size ):
        a = numpy.fromfile( self.file, dtype=dtype, count=size )
        if self.byteswap_needed:
            a.byteswap()
        return a
    
    def read( self, byte_count=1 ):
        return self.file.read( byte_count )
        
    def tell( self ):
        return self.file.tell()
        
    def skip( self, count ):
        self.file.seek( count, 1 )
        
    def seek( self, pos, whence=0 ):
        return self.file.seek( pos, whence )

    def read_uint8( self ):
        return self.read_and_unpack( "B", 1 )[0]
    
    def read_uint16( self ):
        return self.read_and_unpack( "H", 2 )[0]
        
    def read_uint32( self ):
        return self.read_and_unpack( "L", 4 )[0]
        
    def read_uint64( self ):
        return self.read_and_unpack( "Q", 8 )[0]

    def read_float( self ):
        return self.read_and_unpack( "f", 4 )[0]
        
        
class BinaryFileWriter( object ):
    """
    Wrapper for doing binary writes on any file like object.
    
    Currently this is not heavily optimized (it uses the `struct` module to
    unpack)
    """    
    def __init__( self, file, magic = None, is_little_endian = False ):
        self.is_little_endian = is_little_endian
        if self.is_little_endian:
            self.endian_code = "<"
        else:
            self.endian_code = ">"
        self.file = file
        if magic is not None:
            self.write_uint32( magic )
            
    def pack( self, format, buffer ):
        pattern = "%s%s" % ( self.endian_code, format )
        return struct.pack( pattern, buffer )
        
    def pack_and_write( self, format, value ):
        """
        Read enough bytes to unpack according to `format` and return the
        tuple of unpacked values.
        """
        pattern = "%s%s" % ( self.endian_code, format )
        return self.file.write( struct.pack( pattern, value ) )
    
    def write_c_string( self, value ):
        """
        Read a zero terminated (C style) string
        """
        self.file.write( value )
        self.file.write( '\0' )
        
    def write_raw_array( self, value ):
        value.tofile( self.file )
    
    def write( self, value ):
        return self.file.write( value )
        
    def skip( self, count ):
        self.file.seek( count, 1 )
        
    def tell( self ):
        return self.file.tell()
        
    def seek( self, pos, whence=0 ):
        return self.file.seek( pos, whence )

    def write_uint8( self, value ):
        return self.pack_and_write( "B", value )
    
    def write_uint16( self, value ):
        return self.pack_and_write( "H", value )

    def write_uint32( self, value ):
        return self.pack_and_write( "L", value )

    def write_uint64( self, value ):
        return self.pack_and_write( "Q", value )
