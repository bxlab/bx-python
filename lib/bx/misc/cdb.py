from UserDict import DictMixin

from bx.misc.binary_file import BinaryFileReader, BinaryFileWriter
import numpy
import sys

def cdbhash( s ):
    return reduce( lambda h, c: (((h << 5) + h) ^ ord(c)) & 0xffffffffL, s, 5381 )

class FileCDBDict( DictMixin ):
    """
    For accessing a CDB structure on disk. Read only. Currently only supports
    access by key (getitem).
    
    NOTE: The keys method could be implemented by scanning the main table.
    """
    def __init__( self, file, is_little_endian=True ):
        # TODO: Deal with endianess
        self.io = BinaryFileReader( file, is_little_endian=is_little_endian )
        self.header_offset = self.io.tell()
        # Read the whole header (only 2k)
        self.header = []
        for i in range( 256 ):
            self.header.append( ( self.io.read_uint32(), self.io.read_uint32() ) )
    def __getitem__( self, key ):
        hash = cdbhash( key )
        # Find position of subtable using 8 LSBs of hash
        subtable_offset = self.header[ hash % 256 ][0]
        subtable_size = self.header[ hash % 256 ][1]
        if subtable_size == 0:
            raise KeyError
        # Seek into subtable and look for match
        start = ( hash >> 8 )
        for i in range( subtable_size ):
            offset = subtable_offset + ( ( start + i ) % subtable_size ) * 8
            self.io.seek( offset )
            h = self.io.read_uint32()
            p = self.io.read_uint32()
            # Hit an empty bin, no match for key
            if p == 0:
                raise KeyError
            # Hash matches, need to check full key
            if h == hash:
                self.io.seek( p )
                klen = self.io.read_uint32()
                vlen = self.io.read_uint32()
                k = self.io.read( klen )
                if k == key:
                    v = self.io.read( vlen )
                    return v
        else:
            # Visited every slot and no match (should never happen since
            # there are empty slots by contruction)
            raise KeyError
        
    @classmethod
    def to_file( Class, dict, file, is_little_endian=True ):
        """
        For constructing a CDB structure in a file. Able to calculate size on
        disk and write to a file
        """
        io = BinaryFileWriter( file, is_little_endian=is_little_endian )
        start_offset = io.tell()
        # Header is of fixed length
        io.seek( start_offset + ( 8 * 256 ) )
        # For each item, key and value length (written as length prefixed
        # strings). We also calculate the subtables on this pass.
        # NOTE: This requires the key and value be byte strings, support for
        #       dealing with encoding specific value types should be
        #       added to this wrapper
        subtables = [ [] for i in range(256) ]
        for key, value in dict.iteritems():
            pair_offset = io.tell()
            io.write_uint32( len( key ) )
            io.write_uint32( len( value ) )
            io.write( key )
            io.write( value )
            hash = cdbhash( key )
            subtables[ hash % 256 ].append( ( hash, pair_offset ) )
        # Save the offset where the subtables will start
        subtable_offset = io.tell()
        # Write subtables
        for subtable in subtables:
            if len( subtable ) > 0:
                # Construct hashtable to be twice the size of the number
                # of items in the subtable, and built it in memory
                ncells = len( subtable ) * 2
                cells = [ (0,0) for i in range( ncells ) ]
                for hash, pair_offset in subtable:
                    index = ( hash >> 8 ) % ncells
                    while cells[index][1] != 0:
                        index = ( index + 1 ) % ncells
                    # Guaranteed to find a non-empty cell
                    cells[index] = ( hash, pair_offset )
                # Write subtable
                for hash, pair_offset in cells:
                    io.write_uint32( hash )
                    io.write_uint32( pair_offset )
        # Go back and write the header
        end_offset = io.tell()
        io.seek( start_offset )
        index = subtable_offset
        for subtable in subtables:
            io.write_uint32( index )
            io.write_uint32( len( subtable * 2 ) )
            # For each cell in the subtable, a hash and a pointer to a value
            index += ( len( subtable ) * 2 ) * 8
        # Leave fp at end of cdb
        io.seek( end_offset )