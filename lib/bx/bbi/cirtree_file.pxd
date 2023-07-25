from .types cimport *


cdef class CIRTreeFile:
    cdef object file
    cdef object reader
    cdef boolean is_byteswapped
    cdef bits64 root_offset
    cdef bits32 block_size
    cdef bits64 item_count
    cdef bits32 start_chrom_ix
    cdef bits32 start_base
    cdef bits32 end_chrom_ix
    cdef bits32 end_base
    cdef bits64 file_size
    cdef bits32 items_per_slot
