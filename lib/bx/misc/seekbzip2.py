import os
import bisect
from bx_extras.filelike import FileLikeBase
    
from _seekbzip2 import SeekBzip2
    
class SeekableBzip2File( FileLikeBase ):
    
    def __init__( self, filename, table_filename ):
        FileLikeBase.__init__( self )
        self.filename = filename
        self.table_filename = table_filename
        self.init_table()
        self.init_bz2()
        self.pos = 0
        self.dirty = True
        
    def init_bz2( self ):
        self.seek_bz2 = SeekBzip2( self.filename )
        
    def init_table( self ):
        # Position in plaintext file
        self.table_positions = []
        # Position of corresponding block in bz2 file (bits)
        self.table_bz2positions = []
        pos = 0
        for line in open( self.table_filename ):
            fields = line.split()
            # Position of the compressed block in the bz2 file
            bz2_pos = int( fields[0] )
            # Length of the block when uncompressed
            length = int( fields[1] )
            self.table_positions.append( pos )
            self.table_bz2positions.append( bz2_pos )
            pos = pos + length
        self.size = pos
        
    def _read( self, sizehint=-1 ):
        if self.dirty:
            # Our virtual position in the uncompressed data is out of sync
            # FIXME: If we're moving to a later position that is still in 
            # the same block, we could just read and throw out bytes in the
            # compressed stream, less wasteful then backtracking
            chunk, offset = self.get_chunk_and_offset( self.pos )
            # Get the seek position for that chunk and seek to it
            bz2_seek_pos = self.table_bz2positions[chunk] 
            self.seek_bz2.seek( bz2_seek_pos )
            # Consume bytes to move to the correct position
            assert len( self.seek_bz2.read( offset ) ) == offset
            # Update state
            self.dirty = False
        # Attempt to read desired amount
        if sizehint == -1:
            sizehint = 256 * 1024
        val = self.seek_bz2.read( sizehint )
        if val is None:
            # EOF
            self.pos = self.size
        else:
            self.pos = self.pos + len( val )
        return val
        
    def tell( self ):
        # HACK: This could be handled in a base class -- FilelikeBase does
        #       buffered reading, so we need to account for any byte still
        #       in it's buffer when determing current position.
        return self.pos - len( self._FileLikeBase__rbuffer )
            
    def get_chunk_and_offset( self, position ):
        # Find the chunk that position is in using a binary search
        chunk = bisect.bisect( self.table_positions, position ) - 1
        offset = position - self.table_positions[chunk]
        return chunk, offset
        
    def seek( self, offset, whence=0 ):
        # Determine absolute target position
        if whence == 0:
            target_pos = offset
        elif whence == 1:
            target_pos = self.pos + offset
        elif whence == 2:
            target_pos = self.size - offset
        else:
            raise Exception( "Invalid `whence` argument: %r", whence )
        # Check if this is a noop
        if target_pos == self.pos:
            return    
        # Verify it is valid
        assert 0 <= target_pos < self.size, "Attempt to seek outside file"
        # Move the position
        self.pos = target_pos
        # Mark as dirty, the next time a read is done we need to actually
        # move the position in the bzip2 file
        self.dirty = True
        # HACK: This could be handled in a base class -- FilelikeBase does
        #       buffered reading, so we need to reset the buffer since
        #       it is no longer valid
        self._FileLikeBase__rbuffer = ""