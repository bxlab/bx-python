cdef extern from "Python.h":
    char * PyString_AsString( object )
    object PyString_FromStringAndSize( char *, int )
    int _PyString_Resize( object, int ) except -1

import struct, sys

cdef char * NIB_I2C_TABLE 
NIB_I2C_TABLE = "TCAGNXXXtcagnxxx"

def translate_raw_data( data, int start, int length ):
    """
    Data is a block read from the file that needs to be unpacked, dealing
    with end conditions based on start/length.
    """
    cdef int i 
    cdef char * p_rval
    cdef char * p_data
    # Allocate string to write into
    rval = PyString_FromStringAndSize( NULL, length ) 
    # Get char pointer access to strings
    p_rval = PyString_AsString( rval )
    p_data = PyString_AsString( data )
    i = 0
    # Odd start
    if start & 1: 
        p_rval[i] = NIB_I2C_TABLE[ p_data[0] & 0xF ]
        p_data = p_data + 1
        i = 1
    # Two output values for each input value
    while 1:
        p_rval[i] = NIB_I2C_TABLE[ ( p_data[0] >> 4 ) & 0xF ];
        if i == length: break
        i = i + 1
        p_rval[i] = NIB_I2C_TABLE[ ( p_data[0] >> 0 ) & 0xF ];
        if i == length: break
        i = i + 1
        p_data = p_data + 1
    return rval