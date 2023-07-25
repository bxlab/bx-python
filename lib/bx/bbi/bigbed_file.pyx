"""
BigBed file.
"""

import zlib
from io import BytesIO

import numpy

from bx.intervals.io import GenomicInterval
from bx.misc.binary_file import BinaryFileReader

cimport numpy

from .bbi_file cimport (
    BBIFile,
    BlockHandler,
    SummarizedData,
)
from .cirtree_file cimport CIRTreeFile
from .types cimport bits32

DEF big_bed_sig = 0x8789F2EB

cdef inline int range_intersection( int start1, int end1, int start2, int end2 ):
    return min( end1, end2 ) - max( start1, start2 )

cdef class BigBedBlockHandler( BlockHandler ):
    """
    BlockHandler that parses the block into a series of BED records
    """
    cdef bits32 chrom_id
    cdef bits32 start
    cdef bits32 end
    def __init__( self, bits32 chrom_id, bits32 start, bits32 end ):
        BlockHandler.__init__( self )
        self.chrom_id = chrom_id
        self.start = start
        self.end = end
    cdef handle_block( self, bytes block_data, BBIFile bbi_file ):
        cdef object bytes_io
        cdef int length
        cdef bits32 chrom_id, s, e
        cdef str rest
        # Now we parse the block, which should just be a bunch of BED records
        bytes_io = BytesIO( block_data )
        block_reader = BinaryFileReader( bytes_io, is_little_endian=bbi_file.reader.is_little_endian )
        length = len( block_data )
        while bytes_io.tell() < length:
            chrom_id = block_reader.read_uint32()
            s = block_reader.read_uint32()
            e = block_reader.read_uint32()
            rest = block_reader.read_c_string()
            if chrom_id != self.chrom_id:
                continue
            if s < self.end and e > self.start:
                self.handle_interval_value( s, e, rest )

    cdef handle_interval_value( self, bits32 s, bits32 e, str rest ):
        pass

cdef class SummarizingBlockHandler( BigBedBlockHandler ):
    """
    Accumulates intervals into a SummarizedData
    """
    cdef SummarizedData sd
    def __init__( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        BigBedBlockHandler.__init__( self, chrom_id, start, end )
        # What we will load into
        self.sd = SummarizedData( start, end, summary_size )
        for i in range(summary_size):
            self.sd.min_val[i] = +numpy.inf
        for i in range(summary_size):
            self.sd.max_val[i] = -numpy.inf

    cdef handle_interval_value( self, bits32 s, bits32 e, str rest ):
        # FIXME: Does this really obvious thing actually do what we want?
        #        No... sum_data will end up being the coverage, but min/max/etc are wrong
        self.sd.accumulate_interval_value( s, e, 1 )

cdef class IntervalAccumulatingBlockHandler( BigBedBlockHandler ):
    cdef list intervals
    """
    Accumulates intervals into a list of intervals with values
    """
    def __init__( self, bits32 chrom_id, bits32 start, bits32 end ):
        BigBedBlockHandler.__init__( self, chrom_id, start, end )
        self.intervals = []

    cdef handle_interval_value( self, bits32 s, bits32 e, str rest ):
        self.intervals.append( ( s, e, rest ) )

cdef class BigBedFile( BBIFile ): 
    """
    A "big binary indexed" file whose raw data is in BED format.
    """
    def __init__( self, file=None ):
        BBIFile.__init__( self, file, big_bed_sig, "bigbed" )

    cdef _summarize_from_full( self, bits32 chrom_id, bits32 start, bits32 end, int summary_size ):
        """
        Create summary from full data.
        """
        v = SummarizingBlockHandler( chrom_id, start, end, summary_size )
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
        v = IntervalAccumulatingBlockHandler( chrom_id, start, end )
        self.visit_blocks_in_region( chrom_id, start, end, v )
        rval = []
        # FIXME: Not sure the best way to return, will user GenomicInterval for
        # now. 
        for ( s, e, rest ) in v.intervals:
            fields = [ chrom, str( s ), str( e ) ] + rest.split( "\t" )
            rval.append( GenomicInterval( None, fields, 0, 1, 2, 5, "+" ) )
        return rval


