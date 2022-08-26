from cpython.version cimport PY_MAJOR_VERSION


cdef extern from "Python.h":
    char * PyBytes_AsString( object )
    object PyBytes_FromStringAndSize( char *, Py_ssize_t )

cdef extern from "ctype.h":
    int tolower( int )
    
cdef extern from "string.h":
    void * memset( void *, int, size_t )

import struct
import sys
from bisect import bisect


cdef char* valToNt
valToNt = "TCAG"

def read( file, seq, int fragStart, int fragEnd, bint do_mask ):
    """
    Stolen directly from Jim Kent's twoBit.c
    """
    cdef int packedStart, packedEnd, packByteCount
    cdef int pOff, pStart, pEnd
    cdef int midStart, remainder, partCount
    cdef int i, j, s, e
    cdef char * packed
    cdef char * dna
    cdef char * dna_orig
    cdef char partial
    packedStart = (fragStart>>2);
    packedEnd = ((fragEnd+3)>>2);
    packByteCount = packedEnd - packedStart;
    # Empty string in which to write unpacked DNA

    dna_py = PyBytes_FromStringAndSize(NULL, fragEnd - fragStart)
    dna = PyBytes_AsString( dna_py )

    seek_bytes = seq.sequence_offset+packedStart

    # Read it
    file.seek( seek_bytes )

    packed_py = file.read( packByteCount )
    packed = PyBytes_AsString( packed_py )

    # Handle case where everything is in one packed byte 
    if packByteCount == 1:
        pOff = (packedStart<<2)
        pStart = fragStart - pOff
        pEnd = fragEnd - pOff
        partial = packed[0]
        assert pEnd <= 4
        assert pStart >= 0
        for i from pStart <= i < pEnd:
            dna[0] = valToNt[(partial >> (6-i-i)) & 3]
            dna = dna + 1
    else:
        # Handle partial first packed byte.
        midStart = fragStart;
        remainder = ( fragStart&3 )
        if remainder > 0:
            partial = packed[0]; packed = packed + 1
            partCount = 4 - remainder;
            for i from partCount - 1 >= i >= 0:
                dna[i] = valToNt[ partial & 3 ]
                partial = partial >> 2
            midStart = midStart + partCount
            dna = dna + partCount
        # Handle middle bytes.
        remainder = fragEnd&3
        midEnd = fragEnd - remainder
        i = midStart
        while i < midEnd:
            partial = packed[0]
            packed = packed + 1;
            dna[3] = valToNt[partial&3];
            partial = partial >> 2
            dna[2] = valToNt[partial&3];
            partial = partial >> 2
            dna[1] = valToNt[partial&3];
            partial = partial >> 2
            dna[0] = valToNt[partial&3];
            dna = dna + 4;            
            # Increment
            i = i + 4
            ## sys.stderr.write( "!!!< " + dna_py + " >!!!\n" ); sys.stderr.flush()
        # End
        if remainder > 0:
            partial = packed[0];
            partial = partial >> (8-remainder-remainder)
            for i from remainder - 1 >= i >= 0:
                dna[i] = valToNt[partial&3]
                partial = partial >> 2
    # Restore DNA pointer
    dna = PyBytes_AsString( dna_py )
    # N's
    n_block_count = len( seq.n_block_starts )
    if n_block_count > 0:
        start_ix = bisect( seq.n_block_starts, fragStart ) - 1
        if start_ix < 0: start_ix = 0            
        for i from start_ix <= i < n_block_count:
            s = seq.n_block_starts[i];
            e = s + seq.n_block_sizes[i];
            if (s >= fragEnd):
                break
            if (s < fragStart):
               s = fragStart
            if (e > fragEnd):
               e = fragEnd
            if (s < e):
                memset( dna + s - fragStart, c'N', e - s)
    # Mask
    if do_mask:
        m_block_count = len( seq.masked_block_starts )
        if m_block_count > 0:
            start_ix = bisect( seq.masked_block_starts, fragStart ) - 1
            if start_ix < 0: start_ix = 0    
            for i from start_ix <= i < m_block_count:
                s = seq.masked_block_starts[i];
                e = s + seq.masked_block_sizes[i];
                if (s >= fragEnd):
                    break
                if (s < fragStart):
                   s = fragStart
                if (e > fragEnd):
                   e = fragEnd
                if (s < e):
                    for j from s <= j < e:
                        dna[j-fragStart] = tolower( dna[j-fragStart] )
    if PY_MAJOR_VERSION >= 3:
        return dna_py.decode()
    else:
        return dna_py
