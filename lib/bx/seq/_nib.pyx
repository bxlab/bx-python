cdef extern from "Python.h":
    char * PyString_AsString( object )
    object PyString_FromStringAndSize( char *, int )
    int _PyString_Resize( object, int ) except -1

import struct, sys

cdef char * NIB_I2C_TABLE 
cdef char * NIB_I2C_TABLE_FIRST 
cdef char * NIB_I2C_TABLE_SECOND
#NIB_I2C_TABLE        = "TCAGNXXXtcagnxxx"
NIB_I2C_TABLE_FIRST  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGNNNNNNNNNNNNNNNNXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXttttttttttttttttccccccccccccccccaaaaaaaaaaaaaaaaggggggggggggggggnnnnnnnnnnnnnnnnxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
NIB_I2C_TABLE_SECOND = "TCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxxTCAGNXXXtcagnxxx"

def translate_raw_data( data, int start, int length ):
    """
    Data is a block read from the file that needs to be unpacked, dealing
    with end conditions based on start/length.
    """
    cdef int i, j
    cdef char * p_rval
    cdef unsigned char * p_data
    # Allocate string to write into
    rval = PyString_FromStringAndSize( NULL, length ) 
    # Get char pointer access to strings
    p_rval = PyString_AsString( rval )
    p_data = <unsigned char *> PyString_AsString( data )
    i = 0
    # Odd start
    if start & 1: 
        #p_rval[i] = NIB_I2C_TABLE[ p_data[0] & 0xF ]
        p_rval[i] = NIB_I2C_TABLE_SECOND[ p_data[0] ]
        p_data = p_data + 1
        i = 1
    # Two output values for each input value
    for j from 0 <= j < (length-i)/2:
        #p_rval[i]   = NIB_I2C_TABLE[ ( p_data[0] >> 4 ) & 0xF ];
        #p_rval[i+1] = NIB_I2C_TABLE[ ( p_data[0] >> 0 ) & 0xF ];
        p_rval[i]   = NIB_I2C_TABLE_FIRST [ p_data[0] ]
        p_rval[i+1] = NIB_I2C_TABLE_SECOND[ p_data[0] ]
        i = i + 2
        p_data = p_data + 1
    # Odd end
    if i < length: 
        p_rval[i] = NIB_I2C_TABLE_FIRST[ p_data[0] ]
    return rval
