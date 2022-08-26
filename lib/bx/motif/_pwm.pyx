"""
Extensions used by the `pwm` module.
"""

from cpython.version cimport PY_MAJOR_VERSION


cdef extern from "Python.h":
    int PyBytes_AsStringAndSize(object obj, char **buffer, Py_ssize_t* length) except -1

cdef extern from "numpy/arrayobject.h":
    ctypedef int intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags
    # These might be other types in the actual header depending on platform 
    ctypedef int npy_int16
    ctypedef float npy_float32

def score_string( ndarray matrix, ndarray char_to_index, object string, ndarray rval ):
    """
    Score each position in string `string` using the scoring matrix `matrix`.
    Characters in the string are mapped to columns in the matrix by `char_to_index`
    and the score for each position is stored in `rval`.
    
    matrix *must* be a 2d array of type float32
    char_to_index *must* be a 1d array of type int16
    rval *must* be a 1d array of type float32 and the same length as string
    """
    cdef char *buffer
    cdef Py_ssize_t len
    cdef float score
    cdef int i, j
    cdef int matrix_width = matrix.dimensions[0]
    cdef npy_int16 char_index
    # Get input string as character pointer
    if PY_MAJOR_VERSION >= 3:
        bytes_string = string.encode()
    else:
        bytes_string = string
    PyBytes_AsStringAndSize(bytes_string, &buffer, &len )
    # Loop over each position in the string 
    cdef int stop = len - matrix.dimensions[0] + 1
    for i from 0 <= i < stop:
        score = 0.0
        for j from 0 <= j < matrix_width:
            char_index = ( <npy_int16 *> ( char_to_index.data + buffer[i+j] * char_to_index.strides[0] ) )[0]
            if char_index < 0: 
                break
            score += ( <npy_float32*> ( matrix.data + j * matrix.strides[0] + char_index * matrix.strides[1] ) )[0]
        else:
            ( <npy_float32*> ( rval.data + i * rval.strides[0] ) )[0] = score
            
def score_string_with_gaps( ndarray matrix, ndarray char_to_index, object string, ndarray rval ):
    """
    Score each position in string `string` using the scoring matrix `matrix`.
    Characters in the string are mapped to columns in the matrix by `char_to_index`
    and the score for each position is stored in `rval`.

    matrix *must* be a 2d array of type float32
    char_to_index *must* be a 1d array of type int16
    rval *must* be a 1d array of type float32 and the same length as string
    """
    cdef char *buffer
    cdef Py_ssize_t len
    cdef float score
    cdef int i, j, string_pos
    cdef int matrix_width = matrix.dimensions[0]
    cdef npy_int16 char_index
    # Get input string as character pointer
    if PY_MAJOR_VERSION >= 3:
        bytes_string = string.encode()
    else:
        bytes_string = string
    PyBytes_AsStringAndSize(bytes_string, &buffer, &len )
    # Loop over each position in the string 
    cdef int stop = len - matrix.dimensions[0] + 1
    for i from 0 <= i < stop:
        if buffer[i] == '-':
            # Never start scoring at a gap
            continue
        score = 0.0
        string_pos = i
        for j from 0 <= j < matrix_width:
            # Advance to the next non-gap character
            while buffer[string_pos] == '-' and string_pos < len:
                string_pos += 1
            # Ran out of non-gap characters, no more scoring is possible
            if string_pos == len:
                return
            # Find character for position and score
            char_index = ( <npy_int16 *> ( char_to_index.data + buffer[string_pos] * char_to_index.strides[0] ) )[0]
            if char_index < 0: 
                break
            score += ( <npy_float32*> ( matrix.data + j * matrix.strides[0] + char_index * matrix.strides[1] ) )[0]
            # Matched a character, move forward
            string_pos += 1
        else:
            ( <npy_float32*> ( rval.data + i * rval.strides[0] ) )[0] = score
