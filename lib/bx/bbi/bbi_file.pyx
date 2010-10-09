"""
Core implementation for reading UCSC "big binary indexed" files.

There isn't really any specification for the format beyond the code, so this
mirrors Jim Kent's 'bbiRead.c' mostly. 
"""

from bpt_file cimport BPTFile
from cirtree_file cimport CIRTreeFile
from types cimport *

from libc cimport limits

import numpy
cimport numpy

from bx.misc.binary_file import BinaryFileReader
from cStringIO import StringIO
import zlib

# Signatures for bbi related file types

cdef public int big_wig_sig = 0x888FFC26
cdef public int big_bed_sig = 0x8789F2EB

# Some record sizes for parsing
DEF summary_on_disk_size = 32

cdef inline int range_intersection( int start1, int end1, int start2, int end2 ):
    return min( end1, end2 ) - max( start1, start2 )

cdef inline int imax(int a, int b): return a if a >= b else b
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
    cdef bits64 valid_count
    cdef double min_val
    cdef double max_val
    cdef double sum_data
    cdef double sum_squares

cdef class SummarizedData:
    """
    The result of using SummaryBlocks read from the file to produce a 
    aggregation over a particular range and resolution
    """
    cdef public int size
    cdef public numpy.ndarray valid_count
    cdef public numpy.ndarray min_val 
    cdef public numpy.ndarray max_val
    cdef public numpy.ndarray sum_data
    cdef public numpy.ndarray sum_squares
    def __init__( self, int size ):
        self.size = size
        self.valid_count = numpy.zeros( self.size, dtype=numpy.uint64 )
        self.min_val = numpy.zeros( self.size, dtype=numpy.float64 )
        self.max_val = numpy.zeros( self.size, dtype=numpy.float64 )
        self.sum_data = numpy.zeros( self.size, dtype=numpy.float64 )
        self.sum_squares = numpy.zeros( self.size, dtype=numpy.float64 )

cdef class BBIFile:
    """
    A "big binary indexed" file. Stores blocks of raw data and numeric 
    summaries of that data at different levels of aggregation ("zoom levels").
    Generic enough to accommodate both wiggle and bed data. 
    """
    # Probably a PyFileObject, or any seekable file-like
    cdef object file
    # A BinaryFileReader created from file
    cdef object reader
    # The magic number or type signature (whether the file is bigWig or bigBed or...)
    cdef public bits32 magic
    # Is the file byteswapped relative to our native byte order?
    cdef boolean is_byteswapped
    # The index to the chromosomes, an embedded BPT file
    cdef BPTFile chrom_bpt
    # Version number
    cdef public bits16 version
    # Number of zoom levels
    cdef public bits16 zoom_levels
    # Offset to chromosome index
    cdef bits64 chrom_tree_offset
    # Offset to unzoomed data
    cdef bits64 unzoomed_data_offset
    # Offset to unzoomed index
    cdef bits64 unzoomed_index_offset
    # If bed, number of columns
    cdef bits16 field_count
    cdef bits16 defined_field_count
    # Offset to an embedded string containing "AutoSQL" format data that defines the columns
    cdef bits64 as_offset
    # Offset to total summary information (if any)
    cdef bits64 total_summary_offset
    # Size of uncompression buffer, 0 if no compression
    cdef bits32 uncompress_buf_size
    # Unzoomed data index
    cdef CIRTreeFile unzoomed_cir
    # Zoom levels list
    cdef public object level_list

    def __init__( self, file=None, expected_sig=None, type_name=None ):
        if file is not None:
            self.open( file, expected_sig, type_name )

    def open( self, file, expected_sig, type_name ):
        """
        Initialize from an existing bbi file, signature (magic) must be passed
        in since this is generic. 
        """
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

    cpdef _get_chrom_id_and_size( self, char * chrom ):
        """
        Lookup id and size from the chromosome named `chrom`
        """
        bytes = self.chrom_bpt.find( chrom )
        if bytes is not None:
            # The value is two 32 bit uints, use the BPT's reader for checking byteswapping
            chrom_id, chrom_size = self.chrom_bpt.reader.unpack( "II", bytes )
            return chrom_id, chrom_size
        else:
            return None, None

    cpdef _best_zoom_level( self, int desired_reduction ):
        cdef ZoomLevel level, closest_level
        cdef int diff, closest_diff = limits.INT_MAX
        if desired_reduction <= 1:
            return None
        for level in self.level_list:
            diff = desired_reduction - level.reduction_level
            if diff >= 0 and diff < closest_diff:
                closest_diff = diff
                closest_level = level
        return closest_level

    cpdef summarize( self, char * chrom, bits32 start, bits32 end, int summary_size ):
        """
        Create a summary of size `size` over the regions `chrom`:`start`-`end`.
        """
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( chrom )
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

    cdef _summarize_from_full( self, chrom_id, start, end, summary_size ):
        """
        Create summary from full data. This is data specific so must be overridden.
        """
        raise TypeError( "Not implemented" )

cdef class ZoomLevel:
    cdef BBIFile bbi_file
    cdef public bits32 reduction_level
    cdef bits32 reserved
    cdef public bits64 data_offset
    cdef public bits64 index_offset

    def _summary_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end ):
        """
        Return a list of all SummaryBlocks that overlap the region 
        `chrom_id`:`start`-`end`
        """
        cdef CIRTreeFile ctf
        cdef int item_count
        rval = []
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
            block_reader = BinaryFileReader( StringIO( block_data ), is_little_endian=reader.is_little_endian )
            for i from 0 <= i < item_count:
                ## NOTE: Look carefully at bbiRead again to be sure the endian
                ##       conversion here is all correct. It looks like it is 
                ##       just pushing raw data into memory and not swapping
                summary = SummaryBlock()
                summary.chrom_id = block_reader.read_uint32()
                summary.start = block_reader.read_uint32()
                summary.end = block_reader.read_uint32()
                summary.valid_count = block_reader.read_uint32()
                summary.min_val = block_reader.read_float()
                summary.max_val = block_reader.read_float()
                summary.sum_data = block_reader.read_float()
                summary.sum_squares = block_reader.read_float()
                rval.append( summary )
        return rval
            
    cdef _summarize( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Summarize directly from file. 

        Looking at Jim's code, it appears that 
          - bbiSummariesInRegion returns all summaries that span start-end in 
            sorted order
          - bbiSummarySlice is then used to aggregate over the subset of those 
            summaries that overlap a single summary element

        Here we try to do this in one pass without loading the list of summaries 
        into an intermediate form
        """
        cdef CIRTreeFile ctf
        cdef int item_count
        cdef bits32 b_chrom_id, b_start, b_end, b_valid_count
        cdef bits32 s_start, s_end, s_valid_count
        cdef bits32 p_start, p_end
        cdef float b_min_val, b_max_val, b_sum_data, b_sum_squares
        # Factor to map base positions to array indexes, also size in bases of 
        # a single summary value
        cdef float base_to_index_factor = ( end - start ) / <float> summary_size
        cdef int overlap
        cdef double overlap_factor
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
            block_reader = BinaryFileReader( StringIO( block_data ), is_little_endian=reader.is_little_endian )
            for i from 0 <= i < item_count:
                ### NOTE: Look carefully at bbiRead again to be sure the endian
                ###       conversion here is all correct. It looks like it is 
                ###       just pushing raw data into memory and not swapping
                # Read the summary block 
                b_chrom_id = block_reader.read_uint32()
                b_start = block_reader.read_uint32()
                b_end = block_reader.read_uint32()
                b_valid_count = block_reader.read_uint32()
                b_min_val = block_reader.read_float()
                b_max_val = block_reader.read_float()
                b_sum_data = block_reader.read_float()
                b_sum_squares = block_reader.read_float()
                # Stop if the summary block is completely past the region of interest
                if b_start > end:
                    break
                # Determine the range of summary pixels overlapped by this summary, 
                s_start = <int> ( ( b_start - start ) / base_to_index_factor )
                s_end = <int> ( ( ( b_end - start ) / base_to_index_factor ) + 1 )
                # Summary might be larger than our entire array, bound
                s_start = imax( s_start, 0 )
                s_end = imin( s_end, summary_size )
                for j from s_start <= j < s_end:
                    p_start = <int> ( j * base_to_index_factor )
                    p_end = <int> ( ( j + 1 ) * base_to_index_factor )
                    overlap = range_intersection( p_start, p_end, b_start, b_end )
                    if overlap > 0:
                        # Summary is weighted by the amount of overlap apparently. This
                        # is not exact!
                        overlap_factor = <double> overlap / ( b_end - b_start )
                        s_valid_count = <bits32> ( b_valid_count * overlap_factor )
                        if s_valid_count > 0:
                            valid_count[j] += s_valid_count
                            sum_data[j] += b_sum_data * overlap_factor
                            sum_squares[j] += b_sum_squares * overlap_factor
                            if max_val[j] < b_max_val:
                                max_val[j] = b_max_val
                            if min_val[j] < b_min_val:
                                min_val[j] = b_min_val
        return rval 

