"""
Classes for index files that map genomic intervals to values
------------------------------------------------------------

:Author: James Taylor (james@bx.psu.edu)
:Version: $Revision: $

Original author James Taylor;  updates by Bob Harris (rsharris@bx.psu.edu)

An interval index file maps genomic intervals to values.

This implementation writes version 1 file format, and reads versions 0 and 1.

index file format:

   All fields are in big-endian format (most significant byte first).

   All intervals are origin-zero, inclusive start, exclusive end.

   The file begins with an index file header, then is immediately followed
   by an index table.  The index table points to index headers, and index
   headers point to bins.  Index headers and bins are referenced via pointers
   (file offsets), and can be placed more or less anywhere in the file.

   File header:

     offset 0x00: 2C FF 80 0A   magic number
     offset 0x04: 00 00 00 01   version (00 00 00 00 is also supported)
     offset 0x08: 00 00 00 2A   (N) number of index sets
     offset 0x0C:  ...          index table

   Index table:
     The index table is a list of N index headers, packed sequentially and
     sorted by name.  The first begins at offset 0x0C.  Each header describes
     one set of intervals.

     offset:      xx xx xx xx   (L) length of index src name
     offset+4:     ...          index src name (e.g. canFam1.chr1)
     offset+4+L:  xx xx xx xx   offset (in this file) to index data
     offset+8+L:  xx xx xx xx   (B) number of bytes in each value;  for version
                                .. 0, this field is absent, and B is assumed
                                .. to be 4

   Index data:
     The index data for (for one index table) consists of the overall range of
     intervals followed by an array of pointers to bins.  The length of the
     array is 1+binForRange(maxEnd-1,maxEnd), where maxEnd is the maximum
     interval end.

     offset:      xx xx xx xx   minimun interval start
     offset+4:    xx xx xx xx   maximum interval end
     offset+8:    xx xx xx xx   offset (in this file) to bin 0
     offset+12:   xx xx xx xx   number of intervals in bin 0
     offset+16:   xx xx xx xx   offset (in this file) to bin 1
     offset+20:   xx xx xx xx   number of intervals in bin 1
      ... and so on ...

   Bin:
     A bin is an array of (start,end,val), sorted by increasing start (with
     end and val as tiebreakers).  Note that bins may be empty (the number of
     intervals indicated in the index data is zero).  Note that B is determined
     from the appropriate entry in the index table.

     offset:      xx xx xx xx   start for interval 1
     offset+4:    xx xx xx xx   end   for interval 1
     offset+8:     ...          (B bytes) value for interval 1
     offset+8+B:  xx xx xx xx   start for interval 2
     offset+12+B: xx xx xx xx   end   for interval 2
     offset+16+B:  ...          (B bytes) value for interval 2
      ... and so on ...
"""

from bisect import *
from struct import *

__all__ = [ 'Indexes', 'Index' ]

MAGIC=0x2cff800a
VERSION=1

MIN=0
MAX=512*1024*1024

BIN_OFFSETS = [ 512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0 ]
BIN_FIRST_SHIFT = 17
BIN_NEXT_SHIFT = 3

def bin_for_range( start, end):
    """Find the smallest bin that can contain interval (start,end)"""
    start_bin, end_bin = start, end - 1
    start_bin >>= BIN_FIRST_SHIFT
    end_bin >>= BIN_FIRST_SHIFT
    for offset in BIN_OFFSETS:
        if start_bin == end_bin:
            return offset + start_bin
        else:
            start_bin >>= BIN_NEXT_SHIFT
            end_bin >>= BIN_NEXT_SHIFT
    raise "Interval (%d,%d) out of range"

class Indexes:
    """A set of indexes, each identified by a unique name"""

    def __init__( self, filename=None ):
        self.indexes = dict()
        if filename is not None: self.open( filename )

    def add( self, name, start, end, val ):
        if name not in self.indexes:
            self.indexes[name] = Index()
        self.indexes[name].add( start, end, val )

    def get( self, name ):
        if self.indexes[name] is None:
            offset, value_size = self.offsets[name]
            self.indexes[name] = Index( filename=self.filename, offset=offset, value_size=value_size )
        return self.indexes[name]

    def find( self, name, start, end ):
        try: return self.get( name ).find( start, end )
        except: return []

    def open( self, filename ):
        self.filename = filename
        self.offsets = dict()  # (will map key to (offset,value_size))
        f = open( filename )
        magic, version, length = read_packed( f, ">3I" )
        if magic != MAGIC:
            raise "File does not have expected header"
        if version > VERSION:
            warn( "File claims version %d, I don't known anything about versions beyond %d. Attempting to continue", version, VERSION )
        for i in range( length ):
            key_len = read_packed( f, ">I" )
            key = f.read( key_len )
            offset = read_packed( f, ">I" )
            if version == 0:
                value_size = 4
            else:
                value_size = read_packed( f, ">I" )
                assert value_size % 4 == 0, "unsupported value size: %s" % value_size
            self.indexes[ key ] = None
            self.offsets[ key ] = (offset,value_size)
        f.close()

    def write( self, f ):
        keys = self.indexes.keys()
        keys.sort()
        # First determine the size of the header
        base = calcsize( ">3I" )
        for key in keys:
            key = str( key )
            base += calcsize( ">I" )
            base += len( key )
            base += calcsize( ">2I" )
        # Now actually write the header
        write_packed( f, ">3I", MAGIC, VERSION, len( self.indexes ) )
        # And write the index table
        for key in keys:
            key = str( key )
            # Write the string prefixed by its length (pascal!)
            write_packed( f, ">I", len( key ) )
            f.write( key )
            # Write offset
            write_packed( f, ">I", base )
            base += self.indexes[key].bytes_required()
            # Write value size
            write_packed( f, ">I", self.indexes[key].value_size )
        # And finally write each index in order
        for key in keys:
            self.indexes[key].write( f )

class Index:

    def __init__( self, min=MIN, max=MAX, filename=None, offset=0, value_size=None ):
        self._value_size = value_size
        self.max_val = 1   # (1, rather than 0, to force value_size > 0)
        if filename is None:
            self.new( min, max )
        else:
            self.open( filename, offset )

    def get_value_size ( self ):
        if self._value_size != None:
            return self._value_size
        else:
            return round_up_to_4( bytes_of( self.max_val ) )
    value_size = property( fget=get_value_size )

    def new( self, min, max ):
        """Create an empty index for intervals in the range min, max"""
        # Ensure the range will fit given the shifting strategy
        assert MIN <= min <= max <= MAX
        self.min = min
        self.max = max
        # Determine the largest bin we will actually use
        self.bin_count = bin_for_range( max - 1, max ) + 1
        # Create empty bins
        self.bins = [ [] for i in range( self.bin_count ) ]

    def open( self, filename, offset ):
        self.filename = filename
        self.offset = offset
        # Open the file and seek to where we expect our header
        f = open( filename )
        f.seek( offset )
        # Read min/max
        min, max = read_packed( f, ">2I" )
        self.new( min, max )
        # Read bin indexes
        self.bin_offsets = []
        self.bin_sizes = []
        for i in range( self.bin_count ):
            o, s = read_packed( f, ">2I" )
            self.bin_offsets.append( o )
            self.bin_sizes.append( s )
        # Initialize bins to None, indicating that they need to be loaded
        self.bins = [ None for i in range( self.bin_count ) ]

    def add( self, start, end, val ):
        """Add the interval (start,end) with associated value val to the index"""
        insort( self.bins[ bin_for_range( start, end ) ], ( start, end, val ) )
        assert val >= 0
        self.max_val = max(self.max_val,val)

    def find( self, start, end ):
        rval = []
        start_bin = ( max( start, self.min ) ) >> BIN_FIRST_SHIFT
        end_bin = ( min( end, self.max ) - 1 ) >> BIN_FIRST_SHIFT
        for offset in BIN_OFFSETS:
            for i in range( start_bin + offset, end_bin + offset + 1 ):
                if self.bins[i] is None: self.load_bin( i )
                # Iterate over bin and insert any overlapping elements into return value
                for el_start, el_end, val in self.bins[i]:
                    if el_start < end and el_end > start:
                        insort_right( rval, ( el_start, el_end, val ) )
            start_bin >>= BIN_NEXT_SHIFT
            end_bin >>= BIN_NEXT_SHIFT
        return rval

    def iterate( self ):
        for i in range( self.bin_count ):
            if self.bins[i] is None: self.load_bin( i )
            for entry in self.bins[i]:  yield entry

    def load_bin( self, index ):
        bin = []
        if self.bin_sizes[index] == 0:
            self.bins[index] = bin
            return
        f = open( self.filename )
        f.seek( self.bin_offsets[index] )
        # One big read for happy NFS
        item_size = self.value_size + calcsize( ">2I" )
        buffer = f.read( self.bin_sizes[index] * item_size )
        for i in range( self.bin_sizes[index] ):
            start, end = unpack( ">2I", buffer[ i*item_size : i*item_size+8 ] )
            val = unpack_uints( buffer[ i*item_size+8 : (i+1)*item_size ] )
            bin.append( (start, end, val) )
        self.bins[index] = bin
        f.close()

    def write( self, f ):
        value_size = self.value_size
        item_size = value_size + calcsize( ">2I" )
        # Write min/max
        write_packed( f, ">2I", self.min, self.max )
        # Write table of bin sizes and offsets
        base = f.tell() + self.bin_count * calcsize( ">2I" )
        for bin in self.bins:
            write_packed( f, ">2I", base, len( bin ) )
            base += len( bin ) * item_size
        # Write contents of each bin
        for bin in self.bins:
            for start, end, val in bin:
                write_packed( f, ">2I", start, end )
                write_packed_uints( f, val, value_size )

    def bytes_required( self ):
        item_size = self.value_size + calcsize( ">2I" )
        rval = calcsize( ">2I" )
        rval += self.bin_count * calcsize( ">2I" )
        for bin in self.bins:
            rval += len( bin ) * item_size
        return rval

def write_packed( f, pattern, *vals ):
    f.write( pack( pattern, *vals ) )

def read_packed( f, pattern ):
    rval = unpack( pattern, f.read( calcsize( pattern ) ) )
    if len( rval ) == 1: return rval[0]
    return rval

def write_packed_uints( f, v, num_bytes ):
    if num_bytes < 4:
        write_packed( f, ">I", v )
    else:
        parts = []
        while num_bytes > 0:
            parts.append( v & 0xFFFFFFFF )
            v >>= 32
            num_bytes -= 4
        parts.reverse() # (write most-significant chunk first)
        write_packed( f, ">%dI" % len( parts ), *parts )

def unpack_uints( parts ):
    chunks = len( parts )/4
    vals = unpack( ">%dI" % chunks, parts )
    val = vals[0]
    for v in vals[1:]:
        val = (val << 32) + v
    return val

def bytes_of( v ):
    assert v > 0
    b = 0
    while v > 0:
        v >>= 8
        b += 1
    return b

def round_up_to_4( v ):
    if v % 4 == 0:
        return v
    else:
        return v + 4 - (v % 4)

