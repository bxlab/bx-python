cdef enum linemode:
    MODE_BED
    MODE_FIXED
    MODE_VARIABLE

cdef class WiggleReader:
    cdef object file
    cdef object current_chrom
    cdef long current_pos
    cdef long current_step
    cdef long current_span
    cdef linemode mode
  