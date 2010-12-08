"""
BigWig file.
"""

from bbi_file cimport *
from cirtree_file cimport CIRTreeFile
import numpy
cimport numpy
from types cimport *
from bx.misc.binary_file import BinaryFileReader
from cStringIO import StringIO
import zlib, math

DEF big_wig_sig = 0x888FFC26
DEF bwg_bed_graph = 1
DEF bwg_variable_step = 2
DEF bwg_fixed_step = 3

cdef inline int range_intersection( int start1, int end1, int start2, int end2 ):
    return min( end1, end2 ) - max( start1, start2 )

cdef class BigWigFile( BBIFile ): 
    """
    A "big binary indexed" file whose raw data is in wiggle format.
    """
    def __init__( self, file=None ):
        BBIFile.__init__( self, file, big_wig_sig, "bigwig" )
    
    cdef _get_interval_slice( self, bits32 base_start, bits32 base_end, intervals ):
        cdef float valid_count = 0.0
        cdef float sum_data = 0.0
        cdef float sum_squares = 0.0
        cdef float min_val = 0.0
        cdef float max_val = 0.0
        cdef int s, e, overlap
        cdef float val, overlap_factor

        if len(intervals) > 0:
            min_val = intervals[0][2]
            max_val = intervals[0][2]

            for interval in intervals:
                s, e, val = interval

                if s >= base_end:
                    break

                overlap = range_intersection( base_start, base_end, s, e )
                if overlap > 0:
                    interval_size = e - s
                    overlap_factor = <double> overlap / interval_size
                    interval_weight = interval_size * overlap_factor

                    valid_count += interval_weight
                    sum_data += val * interval_weight
                    sum_squares += val * val * interval_weight

                    if max_val < val:
                        max_val = val
                    if min_val > val:
                        min_val = val

        return round(valid_count), sum_data, sum_squares, min_val, max_val

    cdef _summarize_from_full( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Create summary from full data.
        """
        cdef CIRTreeFile ctf
        cdef int item_count
        cdef bits32 b_chrom_id, b_start, b_end, b_valid_count
        cdef bits32 b_item_step, b_item_span
        cdef bits16 b_item_count
        cdef UBYTE b_type
        cdef bits32 base_start, base_end, end1
        cdef int s, e
        cdef float val
        # Factor to map base positions to array indexes, also size in bases of 
        # a single summary value
        cdef float base_to_index_factor = ( end - start ) / <float> summary_size
        # We locally cdef the arrays so all indexing will be at C speeds
        cdef numpy.ndarray[numpy.uint64_t] valid_count
        cdef numpy.ndarray[numpy.float64_t] min_val
        cdef numpy.ndarray[numpy.float64_t] max_val
        cdef numpy.ndarray[numpy.float64_t] sum_data
        cdef numpy.ndarray[numpy.float64_t] sum_squares
        cdef list intervals = []
        # What we will load into
        rval = SummarizedData( summary_size )
        valid_count = rval.valid_count

        min_val = rval.min_val
        for i in range(summary_size):
            min_val[i] = +numpy.inf

        max_val = rval.max_val
        for i in range(summary_size):
            max_val[i] = -numpy.inf

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

                if s < start: 
                    s = start
                if e > end: 
                    e = end
                if s >= e: 
                    continue

                intervals.append([s, e, val])

        base_start = start
        baseCount = end - start

        for i in range(summary_size):
            base_end = start + baseCount*(i+1)/summary_size
            end1 = base_end
            if (end1 == base_start):
                end1 = base_start + 1

            intervals = [interval for interval in intervals if interval[1] > base_start]
            valid_count[i], sum_data[i], sum_squares[i], min_val[i], max_val[i] = self._get_interval_slice(base_start, end1, intervals)
            base_start = base_end

        return rval
