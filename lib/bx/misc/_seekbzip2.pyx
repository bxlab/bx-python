"""
Pyrex/C extension supporting `bx.misc.seekbzip2` (wrapping the low level
functions in `micro-bunzip.c`).
"""

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
        unsigned int * dbuf
    unsigned int get_bits(bunzip_data *bd, char bits_wanted)
    int get_next_block( bunzip_data *bd )
    int read_bunzip(bunzip_data *bd, char *outbuf, int len)
    int start_bunzip(bunzip_data **bdp, int in_fd, char *inbuf, int len)
    int read_bunzip_to_char(bunzip_data *bd, char *outbuf, int len, int* gotcount_out, char stopchar )
    
cdef extern from "unistd.h":
    # Not really
    ctypedef unsigned long long off_t
    off_t lseek( int fildes, off_t offset, int whence )
    
cdef extern from "stdlib.h":
    void free( void *ptr )
 
import sys
import os
    
cdef class SeekBzip2:

    cdef bunzip_data * bd
    cdef int file_fd
    cdef int at_eof

    def __init__( self, filename ):
        self.at_eof = 0
        self.file_fd = os.open( filename, os.O_RDONLY )
        # Initialize bunzip_data from the file
        start_bunzip( &( self.bd ), self.file_fd, NULL, 0 )
        
    def close( self ):
        free( self.bd.dbuf )
        free( self.bd )
        os.close( self.file_fd )

    def seek( self, unsigned long long position ):
        """
        Seek the bunzip_data to a specific chunk (position must correspond to
        that start of a compressed data block).
        """
        cdef off_t n_byte
        cdef int n_bit
        # Break position into bit and byte offsets
        ## sys.stderr.write( "arg pos: %d\n" % position )
        n_byte = position / 8;
        n_bit = position % 8;
        ## sys.stderr.write( "byte pos: %d\n" % n_byte )
        ## sys.stderr.write( "bit pos: %d\n" % n_bit )
        ## sys.stderr.flush()
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
        
    def readline( self, int amount ):
        cdef object rval
        cdef char * p_rval
        cdef int gotcount
        cdef int totalcount
        cdef int status
        cdef int spaceleft
        cdef int desired
        gotcount = 0
        totalcount = 0
        # If already at EOF return None
        if self.at_eof:
            return None
        chunks = []
        # We have great difficulty resizing buffers, so we'll just create
        # one 8k string at a time
        rval = PyString_FromStringAndSize( NULL, 8192 )
        p_rval = PyString_AsString( rval )
        spaceleft = 8192
        while amount != 0:
            if amount > 0 and amount < spaceleft:
                desired = amount
            else:
                desired = spaceleft
            ## sys.stderr.write( "readline, amount: %d\n" % amount )
            ## sys.stderr.write( "buffer: %r" % rval[:100] )
            ## sys.stderr.write( "\n" )
            ## sys.stderr.flush()
            # ord( "\n" ) = 10
            status = read_bunzip_to_char( self.bd, p_rval, desired, &gotcount, 10 );
            ## sys.stderr.write( "readline, desired: %d, gotcount: %d\n" % ( desired, gotcount ) ); 
            ## sys.stderr.write( "buffer: %r" % rval[:100] )
            ## sys.stderr.write( "\n" )
            ## sys.stderr.flush()
            if status == -9: 
                ## sys.stderr.write( "readline, STOP_CHAR\n" ); sys.stderr.flush()
                # Reached the stop character (RETVAL_STOPCHAR == -9), so 
                # we can stop
                chunks.append( rval[:8192-spaceleft+gotcount] )
                break
            elif status == -10:
                ## sys.stderr.write( "readline, BUFFER_FULL\n" ); sys.stderr.flush()
                # Filled the buffer (RETVAL_BUFFER_FULL == -10), so create
                # new buffer and keep going
                chunks.append( rval )
                amount = amount - gotcount
                if amount == 0:
                    # Got the desired amount
                    break
                rval = PyString_FromStringAndSize( NULL, 8192 )
                p_rval = PyString_AsString( rval )
                spaceleft = 8192
            elif status == -8:
                ## sys.stderr.write( "readline, END_OF_BLOCK\n" ); sys.stderr.flush()
                # No more data in the decomp buffer (RETVAL_END_OF_BLOCK == -10)
                if gotcount and p_rval[ gotcount - 1 ] == 10:
                    chunks.append( rval[:8192-spaceleft+gotcount] )
                    break
                # Update buffer info
                p_rval = p_rval + gotcount
                spaceleft = spaceleft - gotcount
                amount = amount - gotcount
                # Get the next block
                status = get_next_block( self.bd )
                if status == -1:
                    # Block is end of stream block (RETVAL_LAST_BLOCK == -1)
                    self.at_eof = 1
                    chunks.append( rval[:gotcount] )
                    break
                self.bd.writeCRC = 0xffffffff
                self.bd.writeCopies = 0
            else:
                # Some other status
                raise Exception( "read_bunzip error %d" % status )
        # Return whatever we read
        return "".join( chunks )
        
    def read( self, int amount ):
        cdef object rval
        cdef char * p_rval
        cdef int gotcount
        cdef int totalcount
        cdef int status
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
    
    