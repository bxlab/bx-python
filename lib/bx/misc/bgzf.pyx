"""
Seekable access to BGZ files based on samtools code. Does not yet implement
complete file-like interface.
"""

ctypedef unsigned long long int64_t

cdef extern from "Python.h":
    char * PyUnicode_AsUTF8( object )
    object PyUnicode_AsUTF8AndSize( char *, int )

cdef extern from "bgzf.h":
    ctypedef struct BGZF
    BGZF * bgzf_open( char * path, char * mode )
    int bgzf_close( BGZF * fp )
    int bgzf_read( BGZF * fp, void * data, int length )
    int64_t bgzf_tell( BGZF * fp )
    int64_t bgzf_seek( BGZF * fp, int64_t pos, int where )
    
cdef class BGZFFile( object ):
    cdef BGZF * bgzf
    def __init__( self, path, mode="r" ):
        self.bgzf = bgzf_open( path, mode )
        if not self.bgzf:
            raise IOError( "Could not open file" )
    def close( self ):
        if self.bgzf:
            bgzf_close( self.bgzf )
    def read( self, int length ):
        cdef object rval = PyUnicode_AsUTF8AndSize( NULL, length )
        bgzf_read( self.bgzf, PyUnicode_AsUTF8( rval ), length )
        return rval
    def tell( self ):
        return bgzf_tell( self.bgzf )
    def seek( self, int64_t pos, int where=0 ):
        return bgzf_seek( self.bgzf, pos, where )