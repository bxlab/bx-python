"""
Pyrex/C extension for quickly finding potential CpG sites in pairs of 
sequences.
"""

from cpython.version cimport PY_MAJOR_VERSION


cdef extern from "find_cpg.h":
    int next_cpg( char * sp1, char * sp2, int start)
    int next_cpg_restricted( char * sp1, char *sp2, int start)
    int next_non_cpg( char * sp1, char * sp2, int start)

def find_cpg( sp1, sp2, start ):
    cdef char* a
    cdef char* b
    cdef int pos
    if PY_MAJOR_VERSION >= 3:
        bytes_sp1, bytes_sp2 = sp1.encode(), sp2.encode()
    else:
        bytes_sp1, bytes_sp2 = sp1, sp2
    a = bytes_sp1
    b = bytes_sp2
    pos = start
    if pos > len(sp1): return -1
    return next_cpg( a, b, pos )

def find_cpg_restricted( sp1, sp2, start ):
    cdef char* a
    cdef char* b
    cdef int pos
    if PY_MAJOR_VERSION >= 3:
        bytes_sp1, bytes_sp2 = sp1.encode(), sp2.encode()
    else:
        bytes_sp1, bytes_sp2 = sp1, sp2
    a = bytes_sp1
    b = bytes_sp2
    pos = start
    if pos > len(sp1): return -1
    return next_cpg_restricted( a, b, pos )

def find_non_cpg( sp1, sp2, start ):
    cdef char* a
    cdef char* b
    cdef int pos
    if PY_MAJOR_VERSION >= 3:
        bytes_sp1, bytes_sp2 = sp1.encode(), sp2.encode()
    else:
        bytes_sp1, bytes_sp2 = sp1, sp2
    a = bytes_sp1
    b = bytes_sp2
    pos = start
    if pos > len(sp1): return -1
    return next_non_cpg( a, b, pos )

def list_cpg( sp1, sp2 ):
    cdef char * a
    cdef char * b
    cdef int start
    if PY_MAJOR_VERSION >= 3:
        bytes_sp1, bytes_sp2 = sp1.encode(), sp2.encode()
    else:
        bytes_sp1, bytes_sp2 = sp1, sp2
    a = bytes_sp1
    b = bytes_sp2
    start = 0
    cpglist = list()
    while start > -1 and start < len(sp1):
        start = next_cpg( a, b, start )
        if start == -1: break
        cpglist.append(start)
        start = start + 1
    return cpglist

def list_cpg_restricted( sp1, sp2 ):
    cdef char * a
    cdef char * b
    cdef int start
    if PY_MAJOR_VERSION >= 3:
        bytes_sp1, bytes_sp2 = sp1.encode(), sp2.encode()
    else:
        bytes_sp1, bytes_sp2 = sp1, sp2
    a = bytes_sp1
    b = bytes_sp2
    start = 0
    cpglist = list()
    while start > -1 and start < len(sp1):
        start = next_cpg_restricted( a, b, start )
        if start == -1: break
        cpglist.append(start)
        start = start + 1
    return cpglist

def list_non_cpg( sp1, sp2 ):
    cdef char * a
    cdef char * b
    cdef int start
    if PY_MAJOR_VERSION >= 3:
        bytes_sp1, bytes_sp2 = sp1.encode(), sp2.encode()
    else:
        bytes_sp1, bytes_sp2 = sp1, sp2
    a = bytes_sp1
    b = bytes_sp2
    start = 0
    cpglist = list()
    while start > -1 and start < len(sp1):
        start = next_non_cpg( a, b, start )
        if start == -1: break
        cpglist.append(start)
        start = start + 1
    return cpglist

def remove_gaps( sp, cpglist ):
    for item in cpglist:
        if sp[item] == '-':
            cpglist.remove(item)
    return cpglist
