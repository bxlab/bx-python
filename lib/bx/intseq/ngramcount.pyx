"""
Tools for counting words (n-grams) in integer sequences.
"""

import numpy

cdef extern from "Python.h":
    ctypedef int Py_intptr_t
    long PyInt_AsLong(object)
    void Py_INCREF(object)
    void Py_DECREF(object)
    void * PyCObject_AsVoidPtr(object cobj)
    
# for PyArrayInterface:
CONTIGUOUS=0x01
FORTRAN=0x02
ALIGNED=0x100
NOTSWAPPED=0x200
WRITEABLE=0x400

ctypedef struct PyArrayInterface:
    int two              # contains the integer 2 as a sanity check
    int nd               # number of dimensions
    char typekind        # kind in array --- character code of typestr
    int itemsize         # size of each element
    int flags            # flags indicating how the data should be interpreted
    Py_intptr_t *shape   # A length-nd array of shape information
    Py_intptr_t *strides # A length-nd array of stride information
    void *data           # A pointer to the first element of the array
    
def count_ngrams( object ints, int n, int radix ):
    """
    Count the number of occurrences of each possible length `n` word in 
    `ints` (which contains values from 0 to `radix`). Returns an array
    of length `radix` ** `n` containing the counts.
    """
    cdef PyArrayInterface * ints_desc, * rval_desc
    # Get array interface for input string and validate
    ints_desc_obj = ints.__array_struct__
    ints_desc = <PyArrayInterface *> PyCObject_AsVoidPtr( ints_desc_obj )
    assert ints_desc.two == 2, "Array interface sanity check failed, got %d" % ints_desc.two
    assert ints_desc.nd == 1, "Input array must be 1d"
    assert ints_desc.typekind == 'i'[0], "Input array must contain integers"
    assert ints_desc.itemsize == 4, "Input array must contain 32bit integers"
    assert ints_desc.flags & CONTIGUOUS > 0, "Input array must be contiguous"
    assert ints_desc.flags & ALIGNED > 0, "Input array must be aligned"
    assert ints_desc.flags & NOTSWAPPED > 0, "Input array must not be byteswapped"
    # Create numpy array for return value, get array interface and validate
    rval = numpy.zeros( <int> ( ( <float> radix ) ** n ), dtype=numpy.int32 )
    assert ints_desc.two == 2, "Array interface sanity check failed, got %d" % ints_desc.two
    rval_desc_obj = rval.__array_struct__
    rval_desc = <PyArrayInterface *> PyCObject_AsVoidPtr( rval_desc_obj )
    assert rval_desc.two == 2, "Array interface sanity check failed"
    assert rval_desc.nd == 1, "Input array must be 1d"
    assert rval_desc.typekind == 'i'[0], "Input array must contain integers"
    assert rval_desc.itemsize == 4, "Input array must contain 32bit integers"
    assert rval_desc.flags & CONTIGUOUS > 0, "Input array must be contiguous"
    assert rval_desc.flags & ALIGNED > 0, "Input array must be aligned"
    assert rval_desc.flags & NOTSWAPPED > 0, "Input array must not be byteswapped"
    # Do it
    _count_ngrams( <int*> ints_desc.data, ints_desc.shape[0], <int*> rval_desc.data, n, radix )
    return rval
    
cdef _count_ngrams( int* ints, int n_ints, int* rval, int n, int radix ):
    cdef int i, j, index, factor, letter
    # Loop over each word in the string
    for i from 0 <= i < ( n_ints - n ):
        # Walk back to build index into count array
        index = 0
        factor = 1
        for j from 0 <= j < n:
            letter = ints[ i + j ]
            if letter < 0 or letter >= radix:
                # This word is bad, break out and do not increment counts
                print "breaking, letter", letter
                break
            index = index + letter * factor
            factor = factor * radix
        else:
            print index
            rval[ index ] = rval[ index ] + 1
            