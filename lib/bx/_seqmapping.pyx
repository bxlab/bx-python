"""
Pyrex extension classes used by `seqmapping.py`.
"""

from cpython.buffer cimport (
    PyBUF_SIMPLE,
    PyBUF_WRITABLE,
    PyBuffer_Release,
    PyObject_GetBuffer,
)
from cpython.bytes cimport PyBytes_AsStringAndSize
from libc.stdlib cimport (
    free,
    malloc,
)

import random
import sys
from math import floor

from numpy import zeros


cdef class CharToIntArrayMapping:
    """Mapping for converting strings to int arrays"""

    cdef int table[256]
    cdef int out_size
    cdef object reverse_table

    def __cinit__( self ):
        """Init empty mapping (all characters map to -1)"""
        cdef int i
        for i from 0 <= i < 256: self.table[i] = -1
        self.out_size = 0

    def __init__( self ):
        self.reverse_table = dict()

    def set_mapping( self, c, int symbol ):
        """Modify mapping so 'chars' map to 'symbol'"""
        char = ord( c )
        self.table[ char ] = symbol
        if self.out_size <= symbol:
            self.out_size = symbol + 1
        self.reverse_table[ symbol ] = chr( char )

    def translate( self, string ):
        """Translate 'string' and return as int array"""
        cdef Py_ssize_t s_len
        cdef unsigned char * s_buf
        cdef Py_buffer t_view
        cdef int * t_buf
        # Get direct access to string
        PyBytes_AsStringAndSize( string, <char **> &s_buf, &s_len )
        # Initialize empty array
        rval = zeros( s_len, 'i' )
        if PyObject_GetBuffer(rval, &t_view, PyBUF_WRITABLE) != 0:
            raise ValueError("Unable to get buffer")
        try:
            t_buf = <int *>t_view.buf
            # Translate
            for i from 0 <= i < s_len:
                t_buf[i] = self.table[s_buf[i]]
        finally:
            PyBuffer_Release(&t_view)
        return rval

    def translate_list( self, strings ):
        """Translate a list of strings into an int array"""
        cdef Py_ssize_t i, s_len, text_len
        cdef Py_buffer t_view
        cdef int factor = 1
        cdef unsigned char * s_buf
        cdef int * t_buf

        # No input, no output
        if len( strings ) < 1: return None

        # Length of result
        text_len = len( strings[0] )

        # Init result array
        rval = zeros( text_len, 'i' )
        if PyObject_GetBuffer(rval, &t_view, PyBUF_WRITABLE) != 0:
            raise ValueError("Unable to get buffer")
        try:
            t_buf = <int *>t_view.buf
            # Loop over seqs and accumulate result values
            for string in strings:
                PyBytes_AsStringAndSize(string, <char **> &s_buf, &s_len)
                for i from 0 <= i < text_len:
                    if t_buf[i] >= 0:
                        if self.table[s_buf[i]] == -1:
                            t_buf[i] = -1
                        else:
                            t_buf[i] += self.table[s_buf[i]] * factor
                factor *= self.out_size
        finally:
            PyBuffer_Release(&t_view)
        return rval

    def reverse_map( self, val, nseqs ):
        factor = self.out_size ** (nseqs-1)
        rval = []
        while factor > 0:
            rval.append( self.reverse_table[ int( floor( val / factor ) ) ] )
            val = val - ( floor(val/factor) * factor )
            factor = floor( factor / self.out_size )
        rval.reverse()
        return rval

    def get_out_size( self ):
        return self.out_size

cdef class IntToIntMapping:

    cdef int* table
    cdef int in_size
    cdef int out_size

    def __cinit__( self, int in_size ):
        self.in_size = in_size
        self.table = <int *> malloc( in_size * sizeof( int ) )
        if self.table == NULL: raise "Malloc Failed"
        for i from 0 <= i < in_size: self.table[i] = -1
        self.out_size = 0

    def __dealloc__( self ):
        # sys.stderr.write( "freeing mapping_helper.IntToIntMapping\n" ); sys.stderr.flush()
        free( self.table )

    def set_mapping( self, int index, int symbol ):
        assert ( -1 <= index < self.in_size ), "%d not between 0 and %s" % ( index, self.in_size )
        self.table[index] = symbol
        if self.out_size <= symbol:
            self.out_size = symbol + 1

    def translate( self, src ):
        """Translate `string` and return as int array"""
        cdef Py_buffer s_view, t_view
        cdef Py_ssize_t s_len
        cdef int *s_buf
        cdef int *t_buf
        # Get direct access to string
        if PyObject_GetBuffer(src, &s_view, PyBUF_SIMPLE) != 0:
            raise ValueError("Unable to get buffer")
        try:
            s_len = s_view.len // sizeof(int)
            s_buf = <int *>s_view.buf
            assert s_len == len(src), "`src` argument must be a buffer of 32bit integers"
            # Initialize empty array
            rval = zeros(s_len, "i")
            if PyObject_GetBuffer(rval, &t_view, PyBUF_WRITABLE) != 0:
                raise ValueError("Unable to get buffer")
            try:
                t_buf = <int *>t_view.buf
                # Translate
                for i from 0 <= i < s_len:
                    if s_buf[i] == -1:
                        t_buf[i] = -1
                    elif s_buf[i] >= self.in_size:
                        t_buf[i] = -1
                    else:
                        t_buf[i] = self.table[s_buf[i]]
            finally:
                PyBuffer_Release(&t_view)
        finally:
            PyBuffer_Release(&s_view)
        return rval

    def __getitem__( self, int x ):
        if x == -1: return -1
        assert 0 <= x < self.in_size
        return self.table[ x ]

    def collapse( self, int a, int b ):
        cdef int i
        cdef IntToIntMapping copy
        copy = IntToIntMapping( self.in_size )
        copy.out_size = self.out_size - 1
        if a > b: a, b = b, a
        for i from 0 <= i < self.in_size:
            if self.table[i] == b: copy.table[i] = a
            elif self.table[i] == copy.out_size: copy.table[i] = b
            else: copy.table[i] = self.table[i]
        return copy

    def expand( self, int x ):
        """Grow the alphabet by making 'a' a seperate symbol. If it already mapped to a single symbol, do nothing"""
        cdef int i, count, a, b
        cdef IntToIntMapping copy
        # Get the symbol x maps to
        a = self.table[x]
        # Symbols that map to -1 should not be touched
        if a < 0: return self
        # Count how many other input symbols map to a
        count = 0
        for i from 0 <= i < self.in_size:
            if self.table[i] == a: count = count + 1
        # Already a singleton
        if count < 2: return self
        # Otherwise, make a copy with the separated symbol
        copy = IntToIntMapping( self.in_size )
        copy.out_size = self.out_size + 1
        for i from 0 <= i < self.in_size:
            copy.table[i] = self.table[i]
        copy.table[x] = self.out_size
        return copy

    def expand_out( self, int a ):
        """Grow the alphabet breaking 'a' into two symbols randomly"""
        cdef int i, count, to_split, b
        cdef IntToIntMapping copy
        count = 0
        for i from 0 <= i < self.in_size:
            if self.table[i] == a: count = count + 1
        if count < 2: return self
        copy = IntToIntMapping( self.in_size )
        copy.out_size = self.out_size + 1
        b = self.out_size
        to_split = random.randrange( count )
        count = 0
        for i from 0 <= i < self.in_size:
            if self.table[i] == a:
                if count == to_split: copy.table[i] = b
                else: copy.table[i] = a
                count = count + 1
            else:
                copy.table[i] = self.table[i]
        return copy

    def expand_random_split( self, int a ):
        """Grow the alphabet breaking 'a' into two symbols randomly"""
        cdef int i, count, b
        cdef IntToIntMapping copy
        count = 0
        for i from 0 <= i < self.in_size:
            if self.table[i] == a: count = count + 1
        if count < 2: return self
        copy = IntToIntMapping( self.in_size )
        copy.out_size = self.out_size + 1
        b = self.out_size
        to_split = random.sample( range( count ), count/2 )
        count = 0
        for i from 0 <= i < self.in_size:
            if self.table[i] == a:
                if count in to_split: copy.table[i] = b
                else: copy.table[i] = a
                count = count + 1
            else:
                copy.table[i] = self.table[i]
        return copy

    def get_in_size( self ):
        return self.in_size

    def get_out_size( self ):
        return self.out_size

    def get_table( self ):
        rval = zeros( self.in_size, 'i' )
        for i in range( self.in_size ):
            rval[i] = self.table[i]
        return rval

