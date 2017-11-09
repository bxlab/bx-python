# cython: profile=False

"""
Core implementation for reading UCSC "big binary indexed" files.

There isn't really any specification for the format beyond the code, so this
mirrors Jim Kent's 'bbiRead.c' mostly. 
"""

from cpython.version cimport PY_MAJOR_VERSION
import sys

cimport cython

from collections import deque
from bpt_file cimport BPTFile
from cirtree_file cimport CIRTreeFile
from types cimport *

from libc cimport limits

import numpy
cimport numpy

from bx.misc.binary_file import BinaryFileReader
from io import BytesIO
import zlib, math

cdef extern from "Python.h":
    char * PyBytes_AsString( object )

# Signatures for bbi related file types

cdef public int big_wig_sig = 0x888FFC26
cdef public int big_bed_sig = 0x8789F2EB

# Some record sizes for parsing
DEF summary_on_disk_size = 32

@cython.profile(False)
cdef inline int range_intersection( int start1, int end1, int start2, int end2 ):
    return min( end1, end2 ) - max( start1, start2 )

@cython.profile(False)
cdef inline int imax(int a, int b): return a if a >= b else b
@cython.profile(False)
cdef inline int imin(int a, int b): return a if a <= b else b

cdef enum summary_type:
    summary_type_mean = 0
    summary_type_max = 1
    summary_type_min = 2
    summary_type_coverage = 3
    summary_type_sd = 4

cdef class SummaryBlock:
    """
    A block of summary data from disk
    """
    pass

cdef class SummarizedData:
    """
    The result of using SummaryBlocks read from the file to produce a 
    aggregation over a particular range and resolution
    """
    def __init__( self, bits32 start, bits32 end, int size ):
        self.start = start
        self.end = end
        self.size = size
        self.valid_count = numpy.zeros( self.size, dtype=numpy.float64 )
        self.min_val = numpy.zeros( self.size, dtype=numpy.float64 )
        self.max_val = numpy.zeros( self.size, dtype=numpy.float64 )
        self.sum_data = numpy.zeros( self.size, dtype=numpy.float64 )
        self.sum_squares = numpy.zeros( self.size, dtype=numpy.float64 )
    cdef accumulate_interval_value( self, bits32 s, bits32 e, float val ):
        cdef int base_start, base_end, base_step, overlap, j, interval_size
        cdef double overlap_factor, interval_weight
        # We locally cdef the arrays so all indexing will be at C speeds
        cdef numpy.ndarray[numpy.float64_t, ndim=1] valid_count = self.valid_count
        cdef numpy.ndarray[numpy.float64_t, ndim=1] min_val = self.min_val
        cdef numpy.ndarray[numpy.float64_t, ndim=1] max_val = self.max_val
        cdef numpy.ndarray[numpy.float64_t, ndim=1] sum_data = self.sum_data
        cdef numpy.ndarray[numpy.float64_t, ndim=1] sum_squares = self.sum_squares
        # Trim interval down to region of interest
        if s < self.start: 
            s = self.start
        if e > self.end: 
            e = self.end
        if s >= e: 
            return
        base_step = ( self.end - self.start ) / self.size
        for j from 0 <= j < self.size:
            base_start = self.start + ( base_step * j )
            base_end = base_start + base_step
            overlap = range_intersection( base_start, base_end, s, e )
            if overlap > 0:
                interval_size = e - s
                overlap_factor = <double> overlap / interval_size
                interval_weight = interval_size * overlap_factor
                valid_count[j] += interval_weight
                sum_data[j] += val * interval_weight
                sum_squares[j] += val * val * interval_weight
                if max_val[j] < val:
                    max_val[j] = val
                if min_val[j] > val:
                    min_val[j] = val 

cdef class BlockHandler:
    """
    Callback for `BBIFile.visit_blocks_in_region`
    """
    cdef handle_block( self, bytes block_data, BBIFile bbi_file ):
        pass

cdef class BBIFile:
    """
    A "big binary indexed" file. Stores blocks of raw data and numeric 
    summaries of that data at different levels of aggregation ("zoom levels").
    Generic enough to accommodate both wiggle and bed data. 
    """

    def __init__( self, file=None, expected_sig=None, type_name=None ):
        if file is not None:
            self.open( file, expected_sig, type_name )

    def open( self, file, expected_sig, type_name ):
        """
        Initialize from an existing bbi file, signature (magic) must be passed
        in since this is generic. 
        """
        assert expected_sig is not None
        self.file = file
        # Open the file in a BinaryFileReader, handles magic and byteswapping
        self.reader = reader = BinaryFileReader( file, expected_sig )
        self.magic = expected_sig
        self.is_byteswapped = self.reader.byteswap_needed
        # Read header stuff
        self.version = reader.read_uint16()
        self.zoom_levels = reader.read_uint16()
        self.chrom_tree_offset = reader.read_uint64()
        self.unzoomed_data_offset = reader.read_uint64()
        self.unzoomed_index_offset = reader.read_uint64()
        self.field_count = reader.read_uint16()
        self.defined_field_count = reader.read_uint16()
        self.as_offset = reader.read_uint64()
        self.total_summary_offset = reader.read_uint64()
        self.uncompress_buf_size = reader.read_uint32()
        # Skip reserved
        reader.seek( 64 )
        # Read zoom headers
        self.level_list = []
        for i from 0 <= i < self.zoom_levels:
            level = ZoomLevel() 
            level.bbi_file = self
            level.reduction_level = reader.read_uint32()
            level.reserved = reader.read_uint32()
            level.data_offset = reader.read_uint64()
            level.index_offset = reader.read_uint64()
            self.level_list.append( level )
        # Initialize and attach embedded BPTFile containing chromosome names and ids
        reader.seek( self.chrom_tree_offset )
        self.chrom_bpt = BPTFile( file=self.file )

    cdef visit_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end, BlockHandler handler ):
        """
        Visit each block from the full data that overlaps a specific region
        """
        cdef CIRTreeFile ctf
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
            handler.handle_block( block_data, self )
        
    cpdef summarize( self, object chrom, bits32 start, bits32 end, int summary_size ):
        """
        Gets `summary_size` data points over the regions `chrom`:`start`-`end`.
        """
        cdef char * cchrom
        if PY_MAJOR_VERSION >= 3:
            bytes_chrom = chrom.encode()
        else:
            bytes_chrom = chrom
        cchrom = PyBytes_AsString(bytes_chrom)
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( cchrom )
        if chrom_id is None:
            return None
        # Return value will be a structured array (rather than an array
        # of summary element structures

        # Find appropriate zoom level
        cdef bits32 base_size = end - start
        cdef int full_reduction = base_size / summary_size
        cdef int zoom = full_reduction / 2
        if zoom < 0:
            zoom = 0
        cdef ZoomLevel zoom_level = self._best_zoom_level( zoom )
        if zoom_level is not None:
            return zoom_level._summarize( chrom_id, start, end, summary_size )
        else:
            return self._summarize_from_full( chrom_id, start, end, summary_size )

    cpdef summarize_from_full( self, char * chrom, bits32 start, bits32 end, int summary_size ):
        """
        Gets `summary_size` data points over the regions `chrom`:`start`-`end`, 
        always using the raw data points
        """
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( chrom )
        if chrom_id is None:
            return None
        # Return value will be a structured array (rather than an array
        # of summary element structures
        return self._summarize_from_full( chrom_id, start, end, summary_size )

    cpdef query( self, object chrom, bits32 start, bits32 end, int summary_size ):
        """
        Provides a different view of summary for region, a list of dictionaries
        with keys: mean, max, min, coverage, std_dev
        """

        if end > 2147483647 or start < 0:
            raise ValueError
        results = self.summarize(chrom, start, end, summary_size)

        if not results:
            return None
        
        rval = []
        for i in range(summary_size):
            sum_data = results.sum_data[i]
            valid_count = results.valid_count[i]
            mean = sum_data / valid_count
            coverage = <double> summary_size / (end - start) * valid_count
            
            # print results.sum_squares[i], sum_data, valid_count
            variance = results.sum_squares[i] - sum_data * sum_data / valid_count
            if valid_count > 1:
                variance /= valid_count - 1
            std_dev = math.sqrt(max(variance, 0))

            rval.append( { "mean": mean, "max": results.max_val[i], "min": results.min_val[i], \
                        "coverage": coverage, "std_dev": std_dev } )
        
        return rval

    cdef _get_chrom_id_and_size( self, char * chrom ):
        """
        Lookup id and size from the chromosome named `chrom`
        """
        bytes = self.chrom_bpt.find( chrom )
        if bytes is not None:
            # The value is two 32 bit uints, use the BPT's reader for checking byteswapping
            assert len( bytes ) == 8
            chrom_id, chrom_size = self.chrom_bpt.reader.unpack( "II", bytes )
            return chrom_id, chrom_size
        else:
            return None, None
    
    cdef _summarize_from_full( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Create summary from full data. This is data specific so must be overridden.
        """
        pass
        
    cdef _best_zoom_level( self, int desired_reduction ):
        if desired_reduction <= 1:
            return None
            
        cdef ZoomLevel level, closest_level
        cdef int diff, closest_diff = limits.INT_MAX
        
        closest_level = None
        for level in self.level_list:
            diff = desired_reduction - level.reduction_level
            if diff >= 0 and diff < closest_diff:
                closest_diff = diff
                closest_level = level
        return closest_level

cdef class ZoomLevel:
    cdef BBIFile bbi_file
    cdef public bits32 reduction_level
    cdef bits32 reserved
    cdef public bits64 data_offset
    cdef public bits64 index_offset
    cdef int item_count

    def _summary_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end ):
        """
        Return a list of all SummaryBlocks that overlap the region 
        `chrom_id`:`start`-`end`
        """
        cdef CIRTreeFile ctf
        cdef SummaryBlock summary
        rval = deque()
        reader = self.bbi_file.reader
        reader.seek( self.index_offset )
        ctf = CIRTreeFile( reader.file )
        block_list = ctf.find_overlapping_blocks( chrom_id, start, end )
        for offset, size in block_list:
            # Seek to and read all data for the block
            reader.seek( offset )
            block_data = reader.read( size )
            # Might need to uncompress
            if self.bbi_file.uncompress_buf_size > 0:
                ## block_data = zlib.decompress( block_data, buf_size = self.bbi_file.uncompress_buf_size )
                block_data = zlib.decompress( block_data )
            block_size = len( block_data )
            # The block should be a bunch of summaries. 
            assert block_size % summary_on_disk_size == 0
            item_count = block_size / summary_on_disk_size
            # Create another reader just for the block, shouldn't be too expensive
            block_reader = BinaryFileReader( BytesIO( block_data ), is_little_endian=reader.is_little_endian )
            for i from 0 <= i < item_count:
                ## NOTE: Look carefully at bbiRead again to be sure the endian
                ##       conversion here is all correct. It looks like it is 
                ##       just pushing raw data into memory and not swapping
                
                sum_chrom_id = block_reader.read_uint32()
                # A block can contain summaries from more that one chrom_id
                if sum_chrom_id != chrom_id:
                    block_reader.skip(7*4)
                    continue
                       
                summary = SummaryBlock()
                summary.chrom_id = sum_chrom_id
                summary.start = block_reader.read_uint32()
                summary.end = block_reader.read_uint32()
                summary.valid_count = block_reader.read_uint32()
                summary.min_val = block_reader.read_float()
                summary.max_val = block_reader.read_float()
                summary.sum_data = block_reader.read_float()
                summary.sum_squares = block_reader.read_float()
                rval.append( summary )
                

        return rval
    
    cdef _get_summary_slice( self, bits32 base_start, bits32 base_end, summaries ):
        cdef float valid_count = 0.0
        cdef float sum_data = 0.0
        cdef float sum_squares = 0.0
        cdef float min_val = numpy.nan
        cdef float max_val = numpy.nan
        cdef float overlap_factor
        cdef int overlap
        
        if summaries:
            
            min_val = summaries[0].min_val
            max_val = summaries[0].max_val
            
            for summary in summaries:
                if summary.start >= base_end:
                    break
            
                overlap = range_intersection( base_start, base_end, summary.start, summary.end )
                if overlap > 0:
                    overlap_factor = <double> overlap / (summary.end - summary.start)
                
                    valid_count += summary.valid_count * overlap_factor
                    sum_data += summary.sum_data * overlap_factor
                    sum_squares += summary.sum_squares * overlap_factor

                    if max_val < summary.max_val:
                        max_val = summary.max_val
                    if min_val > summary.min_val:
                        min_val = summary.min_val

        return valid_count, sum_data, sum_squares, min_val, max_val
        
    cdef _summarize( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Summarize directly from file. 

        Looking at Jim's code, it appears that 
          - bbiSummariesInRegion returns all summaries that span start-end in 
            sorted order
          - bbiSummarySlice is then used to aggregate over the subset of those 
            summaries that overlap a single summary element
        """
        cdef bits32 base_start, base_end, base_step
        
        # We locally cdef the arrays so all indexing will be at C speeds
        cdef numpy.ndarray[numpy.float64_t] valid_count
        cdef numpy.ndarray[numpy.float64_t] min_val
        cdef numpy.ndarray[numpy.float64_t] max_val
        cdef numpy.ndarray[numpy.float64_t] sum_data
        cdef numpy.ndarray[numpy.float64_t] sum_squares
        
        # What we will load into
        rval = SummarizedData( start, end, summary_size )
        valid_count = rval.valid_count
        min_val = rval.min_val
        max_val = rval.max_val
        sum_data = rval.sum_data
        sum_squares = rval.sum_squares
        # First, load up summaries
        reader = self.bbi_file.reader
        reader.seek( self.index_offset )
        summaries = self._summary_blocks_in_region(chrom_id, start, end)

        base_step = (end - start) / summary_size
        base_start = start
        base_end = start
        
        for i in range(summary_size):
            base_end += base_step
            
            while summaries and summaries[0].end <= base_start:
                summaries.popleft()

            valid_count[i], sum_data[i], sum_squares[i], min_val[i], max_val[i] = self._get_summary_slice(base_start, base_end, summaries)
            base_start = base_end
        
        return rval
