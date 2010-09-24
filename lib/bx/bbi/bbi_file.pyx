"""
Core implementation for reading UCSC "big binary indexed" files.

There isn't really any specification for the format beyond the code, so this
mirrors Jim Kent's 'bbiRead.c' mostly. 
"""

from bpt_file cimport BPTFile
from cir_tree cimport CIRTreeFile
from types cimport *

from bx.misc.binary_file import BinaryFileReader
from cStringIO import StringIO
import zlib

# Signatures for bbi related file types

DEF big_wig_sig = 0x888FFC26
DEF big_bed_sig = 0x8789F2EB

# Some record sizes for parsing
DEF summary_on_disk_size = 32

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

    def _get_chrom_id_and_size( self, chrom ):
        """
        Lookup id and size from the chromosome named `chrom`
        """
        bytes = self.chrom_bpt.find( chrom )
        if bytes is not None:
            # The value is two 32 bit uints, use the BPT's reader for checking byteswapping
            chrom_id, chrom_size = self.chrom_bpt.reader.unpack( "II", bytes )
            return chrom_id, chrom_size
        else:
            return None

cdef class ZoomLevel:
    cdef BBIFile bbi_file
    cdef public bits32 reduction_level
    cdef bits32 reserved
    cdef public bits64 data_offset
    cdef public bits64 index_offset

    def summaries_in_region( self, bits32 chrom_id, bits32 start, bits32 end ):
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
                summary = Summary()
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
            


    
cdef class ChromIDSize:
    cdef bits32 chrom_id
    cdef bits32 chrom_size

cdef class ChromInfo:
    cdef object name
    cdef bits32 id
    cdef bits32 size

# Skipping 'struct bbiChromUsage'

cdef enum summary_type:
    mean = 0
    max = 1
    min = 2
    coverage = 3
    sd = 4
    
cdef class Summary:
    cdef public bits32 chrom_id
    cdef public bits32 start
    cdef public bits32 end
    cdef public bits32 valid_count
    cdef public float min_val
    cdef public float max_val
    cdef public float sum_data
    cdef public float sum_squares
    
cdef class Interval:
    cdef bits32 start
    cdef bits32 end
    cdef double val

cdef class SummaryElement:
    cdef bits64 valid_count
    cdef double min_val
    cdef double max_val
    cdef double sum_data
    cdef double sum_squares


