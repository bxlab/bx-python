"""
Seekable access to BGZ files based on samtools code. Does not yet implement
complete file-like interface.
"""

from cpython.version cimport PY_MAJOR_VERSION

ctypedef unsigned long long int64_t

cdef extern from "Python.h":
    char * PyBytes_AsString( object )
    object PyBytes_FromStringAndSize( char *, Py_ssize_t )

cdef extern from "bgzf.h":
    ctypedef struct BGZF
    BGZF * bgzf_open( const char * path, const char * mode )
    int bgzf_close( BGZF * fp )
    int bgzf_read( BGZF * fp, void * data, int length )
    int64_t bgzf_tell( BGZF * fp )
    int64_t bgzf_seek( BGZF * fp, int64_t pos, int where )
    
cdef class BGZFFile( object ):
    cdef BGZF * bgzf
    def __init__( self, path, mode="r" ):
        if PY_MAJOR_VERSION >= 3:
            bytes_path, bytes_mode = path.encode(), mode.encode()
        else:
            bytes_path, bytes_mode = path, mode
        self.bgzf = bgzf_open( bytes_path, bytes_mode )
        if not self.bgzf:
            raise IOError( "Could not open file" )
      
    def close( self ):
        if self.bgzf:
            bgzf_close( self.bgzf )
    def read( self, int length ):
        cdef object rval
        rval = PyBytes_FromStringAndSize( NULL, length )
        bgzf_read( self.bgzf, PyBytes_AsString( rval ), length )
        return rval
    def tell( self ):
        return bgzf_tell( self.bgzf )
    def seek( self, int64_t pos, int where=0 ):
        return bgzf_seek( self.bgzf, pos, where )
