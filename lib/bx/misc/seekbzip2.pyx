cdef extern from "Python.h":
    char * PyString_AsString( object )
    object PyString_FromStringAndSize( char *, int )

cdef extern from "micro-bunzip.h":
    ctypedef struct bunzip_data:
        int in_fd
        int inbufBitCount
        int inbufPos
        int inbufCount
        int writeCount
        unsigned int writeCRC
        int writeCurrent
        int writeCopies
    unsigned int get_bits(bunzip_data *bd, char bits_wanted)
    int get_next_block( bunzip_data *bd )
    int read_bunzip(bunzip_data *bd, char *outbuf, int len)
    int start_bunzip(bunzip_data **bdp, int in_fd, char *inbuf, int len)
    
cdef extern from "unistd.h":
    # Not really
    ctypedef int off_t
    off_t lseek( int fildes, off_t offset, int whence )
 
import sys
import os
import bisect
from bx_extras.filelike import FileLikeBase
    
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
        self.seek_bz2 = _SeekBzip2( self.filename )
        
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
        
    
cdef class _SeekBzip2:

    cdef bunzip_data * bd
    cdef int file_fd
    cdef int at_eof

    def __init__( self, filename ):
        self.at_eof = 0
        self.file_fd = os.open( filename, os.O_RDONLY )
        # Initialize bunzip_data from the file
        start_bunzip( &( self.bd ), self.file_fd, NULL, 0 )

    def seek( self, unsigned long position ):
        """
        Seek the bunzip_data to a specific chunk (position must correspond to
        that start of a compressed data block).
        """
        cdef off_t n_byte
        cdef int n_bit
        # Break position into bit and byte offsets
        n_byte = position / 8;
        n_bit = position % 8;
        # Seek the underlying file descriptor
        if ( lseek( self.file_fd, n_byte, 0 ) != n_byte ):
            raise Exception( "lseek of underlying file failed" )
        # Init the buffer at the right bit position
        self.bd.inbufBitCount = self.bd.inbufPos = self.bd.inbufCount = 0
        get_bits( self.bd, n_bit )
        # This ensures that the next read call will return 0, causing the
        # buffer to be re-initialized
        self.bd.writeCount = -1
        # Reset EOF tracking
        self.at_eof = 0

    def read( self, int amount ):
        cdef object rval
        cdef char * p_rval
        cdef int gotcount
        cdef int totalcount
        cdef int status
        cdef int desired
        cdef int READ_BUF_SIZE
        READ_BUF_SIZE = 8 * 1024
        totalcount = 0
        # If already at EOF return None
        if self.at_eof:
            return None
        # Create a new python string large enough to hold the result
        rval = PyString_FromStringAndSize( NULL, amount )
        p_rval = PyString_AsString( rval )
        # Read into it
        ## sys.stderr.write( "read called, bd.current: %x\n" % self.bd.writeCurrent ); sys.stderr.flush()
        while amount > 0:
            if amount > READ_BUF_SIZE:
                desired = READ_BUF_SIZE
            else:
                desired = amount
            gotcount = read_bunzip( self.bd, p_rval, amount );
            if gotcount < 0:
                raise Exception( "read_bunzip error %d" % gotcount )
            elif gotcount == 0:
                status = get_next_block( self.bd )
                if status == -1:
                    self.at_eof = 1
                    break
                self.bd.writeCRC = 0xffffffff
                self.bd.writeCopies = 0
            else:
                totalcount = totalcount + gotcount
                amount = amount - gotcount
                p_rval = p_rval + gotcount
        # Return whatever we read
        return rval[:totalcount]
    
    