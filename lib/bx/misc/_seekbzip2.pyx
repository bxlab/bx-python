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
    
cdef class SeekBzip2:

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
    
    