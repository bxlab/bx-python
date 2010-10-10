"""
BigWig file.
"""

from bpt_file cimport BPTFile
from cirtree_file cimport CIRTreeFile
from bbi_file cimport BBIFile, SummarizedData
from types cimport *

from libc cimport limits

import numpy
cimport numpy

from bx.misc.binary_file import BinaryFileReader
from cStringIO import StringIO
import zlib

DEF big_wig_sig = 0x888FFC26

DEF bwg_bed_graph = 1
DEF bwg_variable_step = 2
DEF bwg_fixed_step = 3

cdef inline int range_intersection( int start1, int end1, int start2, int end2 ):
    return min( end1, end2 ) - max( start1, start2 )

cdef inline int imax(int a, int b): return a if a >= b else b
cdef inline int imin(int a, int b): return a if a <= b else b

cdef class BigWigFile( BBIFile ): 
    """
    A "big binary indexed" file whose raw data is in wiggle format.
    """
    def __init__( self, file=None ):
        BBIFile.__init__( self, file, big_wig_sig, "bigwig" )

    cdef _summarize_from_full( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Create summary from full data. This is data specific so must be overridden.
        """
        cdef CIRTreeFile ctf
        cdef int item_count
        cdef bits32 b_chrom_id, b_start, b_end, b_valid_count
        cdef bits32 b_item_step, b_item_span
        cdef bits16 b_item_count
        cdef UBYTE b_type
        cdef bits32 s_start, s_end, s_valid_count
        cdef bits32 p_start, p_end
        cdef int s, e
        cdef float val
        # Factor to map base positions to array indexes, also size in bases of 
        # a single summary value
        cdef float base_to_index_factor = ( end - start ) / <float> summary_size
        cdef int overlap
        cdef double overlap_factor
        cdef int interval_size
        cdef double interval_weight
        # We locally cdef the arrays so all indexing will be at C speeds
        cdef numpy.ndarray[numpy.uint64_t] valid_count
        cdef numpy.ndarray[numpy.float64_t] min_val
        cdef numpy.ndarray[numpy.float64_t] max_val
        cdef numpy.ndarray[numpy.float64_t] sum_data
        cdef numpy.ndarray[numpy.float64_t] sum_squares
        # What we will load into
        rval = SummarizedData( summary_size )
        valid_count = rval.valid_count
        min_val = rval.min_val
        max_val = rval.max_val
        sum_data = rval.sum_data
        sum_squares = rval.sum_squares
        # First, load up summaries
        reader = self.reader
        reader.seek( self.unzoomed_index_offset )
        ctf = CIRTreeFile( reader.file )
        block_list = ctf.find_overlapping_blocks( chrom_id, start, end )
        for offset, size in block_list:
            # Seek to and read all data for the block
            reader.seek( offset )
            block_data = reader.read( size )
            # Might need to uncompress
            if self.uncompress_buf_size > 0:
                block_data = zlib.decompress( block_data )
            block_size = len( block_data )
            # Now we parse the block, first the header
            block_reader = BinaryFileReader( StringIO( block_data ), is_little_endian=reader.is_little_endian )
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
                elif b_type == bwg_variable_step:
                    s = b_start + i
                    e = s + b_item_span
                    val = block_reader.read_float()
                else:
                    raise Exception( "Unknown record type" )
                # Clip to region of interest
                if s < start: 
                    s = start
                if e > end: 
                    e = end
                if s >= e: 
                    continue
                # Determine the range of summary pixels overlapped by this value 
                s_start = <int> ( ( s - start ) / base_to_index_factor )
                s_end = <int> ( ( ( e - start ) / base_to_index_factor ) + 1 )
                # Summary might be larger than our entire array, bound
                s_start = imax( s_start, 0 )
                s_end = imin( s_end, summary_size )
                for j from s_start <= j < s_end:
                    p_start = <int> ( j * base_to_index_factor )
                    p_end = <int> ( ( j + 1 ) * base_to_index_factor )
                    overlap = range_intersection( p_start, p_end, s, e )
                    if overlap > 0:
                        interval_size = e - s
                        overlap_factor = <double> overlap / interval_size
                        interval_weight = interval_size / overlap_factor
                        valid_count[j] += interval_weight
                        sum_data[j] += val * interval_weight
                        sum_squares[j] += val * val * interval_weight
                        if max_val[j] < val:
                            max_val[j] = val
                        if min_val[j] > val:
                            min_val[j] = val
        return rval 

