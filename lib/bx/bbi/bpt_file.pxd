from bx.misc.binary_file import BinaryFileReader

from .types cimport *


cdef class BPTFile:
    """
    On disk B+ tree compatible with Jim Kent's bPlusTree.c
    """
    cdef object file
    cdef object reader
    cdef boolean is_byteswapped
    cdef bits32 block_size
    cdef bits32 key_size
    cdef bits32 value_size
    cdef bits64 item_count
    cdef bits64 root_offset
