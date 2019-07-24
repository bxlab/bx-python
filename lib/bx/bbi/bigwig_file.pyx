"""
BigWig file.
"""

from collections import deque
from bbi_file cimport *
from cirtree_file cimport CIRTreeFile
import numpy
cimport numpy
from types cimport *
from bx.misc.binary_file import BinaryFileReader
from io import BytesIO
import zlib

DEF big_wig_sig = 0x888FFC26
DEF bwg_bed_graph = 1
DEF bwg_variable_step = 2
DEF bwg_fixed_step = 3

cdef inline int range_intersection( int start1, int end1, int start2, int end2 ):
    return min( end1, end2 ) - max( start1, start2 )

cdef class BigWigBlockHandler( BlockHandler ):
    """
    BlockHandler that parses the block into a series of wiggle records, and calls `handle_interval_value` for each.
    """
    cdef bits32 start
    cdef bits32 end
    def __init__( self, bits32 start, bits32 end ):
        BlockHandler.__init__( self )
        self.start = start
        self.end = end
    cdef handle_block( self, bytes block_data, BBIFile bbi_file ):
        cdef bits32 b_chrom_id, b_start, b_end, b_valid_count
        cdef bits32 b_item_step, b_item_span
        cdef bits16 b_item_count
        cdef UBYTE b_type
        cdef int s, e
        cdef float val
        # Now we parse the block, first the header
        block_reader = BinaryFileReader( BytesIO( block_data ), is_little_endian=bbi_file.reader.is_little_endian )
        b_chrom_id = block_reader.read_uint32()
        b_start = block_reader.read_uint32()
        b_end = block_reader.read_uint32()
        b_item_step = block_reader.read_uint32()
        b_item_span = block_reader.read_uint32()
        b_type = block_reader.read_uint8()
        block_reader.skip(1)
        b_item_count = block_reader.read_uint16()
        for i from 0 <= i < b_item_count:
            # Depending on the type, s and e are either read or 
            # generate using header, val is always read
            if b_type == bwg_bed_graph: 
                s = block_reader.read_uint32()
                e = block_reader.read_uint32()
                val = block_reader.read_float()
            elif b_type == bwg_variable_step:
                s = block_reader.read_uint32()
                e = s + b_item_span
                val = block_reader.read_float()
            elif b_type == bwg_fixed_step:
                s = b_start + ( i * b_item_span )
                e = s + b_item_span
                val = block_reader.read_float()
            else:
                # FIXME: raise exception???
                #        s, e, val are uninitialized/not updated at this point!
                pass
            if s < self.start: 
                s = self.start
            if e > self.end: 
                e = self.end
            if s >= e: 
                continue
            self.handle_interval_value( s, e, val )

    cdef handle_interval_value( self, bits32 s, bits32 e, float val ):
        pass

cdef class SummarizingBlockHandler( BigWigBlockHandler ):
    """
    Accumulates intervals into a SummarizedData
    """
    cdef SummarizedData sd
    def __init__( self, bits32 start, bits32 end, int summary_size ):
        BigWigBlockHandler.__init__( self, start, end )
        # What we will load into
        self.sd = SummarizedData( start, end, summary_size )
        for i in range(summary_size):
            self.sd.min_val[i] = +numpy.inf
        for i in range(summary_size):
            self.sd.max_val[i] = -numpy.inf

    cdef handle_interval_value( self, bits32 s, bits32 e, float val ):
        self.sd.accumulate_interval_value( s, e, val )

cdef class IntervalAccumulatingBlockHandler( BigWigBlockHandler ):
    cdef list intervals
    """
    Accumulates intervals into a list of intervals with values
    """
    def __init__( self, bits32 start, bits32 end ):
        BigWigBlockHandler.__init__( self, start, end )
        self.intervals = []

    cdef handle_interval_value( self, bits32 s, bits32 e, float val ):
        self.intervals.append( ( s, e, val ) )

cdef class ArrayAccumulatingBlockHandler( BigWigBlockHandler ):
    """
    Accumulates intervals into a list of intervals with values
    """
    cdef numpy.ndarray array
    def __init__( self, bits32 start, bits32 end ):
        BigWigBlockHandler.__init__( self, start, end )
        self.array = numpy.zeros( end - start, dtype=numpy.float32 )
        self.array[...] = numpy.nan

    cdef handle_interval_value( self, bits32 s, bits32 e, float val ):
        cdef numpy.ndarray[ numpy.float32_t, ndim=1 ] array = self.array
        cdef int i
        # Slicing is not optimized by Cython
        for i from s - self.start <= i < e - self.start:
            array[ i ] = val

cdef class BigWigHeaderBlockHandler( BigWigBlockHandler ):
    "Reads and returns headers"
    cdef list headers

    def __init__( self, bits32 start, bits32 end ):
        BigWigBlockHandler.__init__( self, start, end )
        self.headers = []

    cdef handle_block( self, bytes block_data, BBIFile bbi_file ):
        cdef bits32 b_chrom_id, b_start, b_end, b_valid_count
        cdef bits32 b_item_step, b_item_span
        cdef bits16 b_item_count
        cdef UBYTE b_type
        cdef int s, e
        cdef float val
        # parse the block header
        block_reader = BinaryFileReader( BytesIO( block_data ), is_little_endian=bbi_file.reader.is_little_endian )
        b_chrom_id = block_reader.read_uint32()
        b_start = block_reader.read_uint32()
        b_end = block_reader.read_uint32()
        b_item_step = block_reader.read_uint32()
        b_item_span = block_reader.read_uint32()
        b_type = block_reader.read_uint8()
        block_reader.skip(1)
        b_item_count = block_reader.read_uint16()
        self.handle_header( b_start, b_end,  b_item_step, b_item_span, b_type, b_item_count )

    cdef handle_header( self, bits32 start, bits32 end, bits32 step, bits32 span, bits8 type, bits16 itemCount ):
        self.headers.append( ( start, end, step, span, type, itemCount ) )

cdef class BigWigFile( BBIFile ):
    """
    A "big binary indexed" file whose raw data is in wiggle format.
    """
    def __init__( self, file=None ):
        BBIFile.__init__( self, file, big_wig_sig, "bigwig" )

    cdef _summarize_from_full( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Create summary from full data.
        """
        v = SummarizingBlockHandler( start, end, summary_size )
        self.visit_blocks_in_region( chrom_id, start, end, v )
        # Round valid count, in place
        for i from 0 <= i < summary_size:
            v.sd.valid_count[i] = round( v.sd.valid_count[i] )
        return v.sd

    cpdef get( self, char * chrom, bits32 start, bits32 end ):
        """
        Gets all data points over the regions `chrom`:`start`-`end`.
        """
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( chrom )
        if chrom_id is None:
            return None
        v = IntervalAccumulatingBlockHandler( start, end )
        self.visit_blocks_in_region( chrom_id, start, end, v )
        return v.intervals

    cpdef get_as_array( self, char * chrom, bits32 start, bits32 end ):
        """
        Gets all data points over the regions `chrom`:`start`-`end`.
        """
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( chrom )
        if chrom_id is None:
            return None
        v = ArrayAccumulatingBlockHandler( start, end )
        self.visit_blocks_in_region( chrom_id, start, end, v )
        return v.array

    cpdef get_headers( self, char * chrom, bits32 start, bits32 end ):
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( chrom )
        if chrom_id is None:
            return None
        v = BigWigHeaderBlockHandler( start, end )
        self.visit_blocks_in_region( chrom_id, start, end, v )
        return v.headers






