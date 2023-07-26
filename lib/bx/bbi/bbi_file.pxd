import numpy

cimport numpy

from .bpt_file cimport BPTFile
from .cirtree_file cimport CIRTreeFile
from .types cimport (
    bits16,
    bits32,
    bits64,
    boolean,
)


cdef class SummaryBlock:
    """
    A block of summary data from disk
    """
    cdef public bits32 chrom_id
    cdef public bits32 start
    cdef public bits32 end
    cdef public bits32 valid_count
    cdef public double min_val
    cdef public double max_val
    cdef public double sum_data
    cdef public double sum_squares

cdef class SummarizedData:
    """
    The result of using SummaryBlocks read from the file to produce a 
    aggregation over a particular range and resolution
    """
    cdef public bits32 start
    cdef public bits32 end
    cdef public int size
    cdef public numpy.ndarray valid_count
    cdef public numpy.ndarray min_val 
    cdef public numpy.ndarray max_val
    cdef public numpy.ndarray sum_data
    cdef public numpy.ndarray sum_squares

    cdef accumulate_interval_value( self, bits32 s, bits32 e, float val )

cdef class BBIFile

cdef class BlockHandler:
    """
    Callback for `BBIFile.visit_blocks_in_region`
    """
    cdef handle_block( self, bytes block_data, BBIFile bbi_file )

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
    # Zoom levels list
    cdef public object level_list

    cdef visit_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end, BlockHandler handler )
    cdef _get_chrom_id_and_size( self, char * chrom )
    cdef _best_zoom_level( self, int desired_reduction )
    cpdef summarize( self, object chrom, bits32 start, bits32 end, int summary_size )
    cpdef summarize_from_full( self, char * chrom, bits32 start, bits32 end, int summary_size )
    cpdef query( self, object chrom, bits32 start, bits32 end, int summary_size )
    cdef _summarize_from_full( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size )    


