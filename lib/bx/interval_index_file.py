from bisect import *
from struct import *

__all__ = [ 'Indexes', 'Index' ]

MAGIC=0x2cff800a
VERSION=0

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
            self.indexes[name] = Index( filename=self.filename, offset=self.offsets[name] )
        return self.indexes[name]

    def find( self, name, start, end ):
        try: return self.get( name ).find( start, end )
        except: return []

    def open( self, filename ):
        self.filename = filename
        self.offsets = dict()
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
            self.indexes[ key ] = None
            self.offsets[ key ] = offset
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
            base += calcsize( ">I" )
        # Now actually write the header
        write_packed( f, ">3I", MAGIC, VERSION, len( self.indexes ) )
        for key in keys:
            key = str( key )
            # Write the string prefixed by its length (pascal!)
            write_packed( f, ">I", len( key ) )
            f.write( key )
            # Write offset 
            write_packed( f, ">I", base )
            base += self.indexes[key].bytes_required()
        # And finally write each index in order
        for key in keys:
            self.indexes[key].write( f )

class Index:

    def __init__( self, min=MIN, max=MAX, filename=None, offset=0 ):
        if filename is None: 
            self.new( min, max )
        else: 
            self.open( filename, offset )

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
        f = open( self.filename )
        f.seek( self.bin_offsets[index] )
        # One big read for happy NFS
        item_size = calcsize( ">3I" )
        buffer = f.read( self.bin_sizes[index] * item_size )
        for i in range( self.bin_sizes[index] ):
            bin.append( unpack( ">3I", buffer[ i*item_size : (i+1)*item_size ] ) )
        self.bins[index] = bin
        f.close()

    def write( self, f ):
        # Write min/max
        write_packed( f, ">2I", self.min, self.max )
        # Write table of bin sizes and offsets
        base = f.tell() + self.bin_count * calcsize( ">2I" )
        for bin in self.bins:
            write_packed( f, ">2I", base, len( bin ) )
            base += len( bin ) * calcsize( ">3I" )
        # Write contents of each bin
        for bin in self.bins:
            for start, end, val in bin:
                write_packed( f, ">3I", start, end, val )

    def bytes_required( self ):
        rval = calcsize( ">2I" )
        rval += self.bin_count * calcsize( ">2I" )
        for bin in self.bins:
            rval += len( bin ) * calcsize( ">3I" )
        return rval 

def write_packed( f, pattern, *vals ):
    f.write( pack( pattern, *vals ) )

def read_packed( f, pattern ):
    rval = unpack( pattern, f.read( calcsize( pattern ) ) )
    if len( rval ) == 1: return rval[0]
    return rval


