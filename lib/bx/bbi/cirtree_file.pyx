from bx.misc.binary_file import BinaryFileReader

DEF cir_tree_sig = 0x2468ACE0

cdef int ovcmp( bits32 a_hi, bits32 a_lo, bits32 b_hi, bits32 b_lo ):
    if a_hi < b_hi: 
        return 1
    elif a_hi > b_hi:
        return -1
    else:
        if a_lo < b_lo:
            return 1
        elif a_lo > b_lo:
            return -1
        else:
            return 0

cdef overlaps( qchrom, qstart, qend, rstartchrom, rstartbase, rendchrom, rendbase ):
    return ( ovcmp( qchrom, qstart, rendchrom, rendbase ) > 0 ) and \
           ( ovcmp( qchrom, qend, rstartchrom, rstartbase ) < 0 )

cdef class CIRTreeFile:

    def __init__( self, file=None ):
        if file is not None:
            self.attach( file )

    def attach( self, file ):
        """
        Attach to an open file
        """
        self.file = file
        self.reader = reader = BinaryFileReader( file, cir_tree_sig )
        self.is_byteswapped = self.reader.byteswap_needed
        # Header
        self.block_size = reader.read_uint32()
        self.item_count = reader.read_uint64()
        self.start_chrom_ix  = reader.read_uint32()
        self.start_base = reader.read_uint32()
        self.end_chrom_ix = reader.read_uint32()
        self.end_base = reader.read_uint32()
        self.file_size = reader.read_uint64()
        self.items_per_slot = reader.read_uint32()
        # Skip reserved
        reader.read_uint32()
        # Save root
        self.root_offset = reader.tell()

    def r_find_overlapping( self, int level, bits64 index_file_offset, bits32 chrom_ix, bits32 start, bits32 end, object rval, object reader ):
        cdef UBYTE is_leaf
        cdef bits16 child_count
        reader.seek( index_file_offset )
        # Block header
        is_leaf = reader.read_uint8()
        assert is_leaf == 0 or is_leaf == 1
        reader.read_uint8()
        child_count = reader.read_uint16()
        # Read block
        if is_leaf:
            self.r_find_overlapping_leaf( level, chrom_ix, start, end, rval, child_count, reader )
        else:
            self.r_find_overlapping_parent( level, chrom_ix, start, end, rval, child_count, reader )

    def r_find_overlapping_leaf( self, int level, bits32 chrom_ix, bits32 start, bits32 end, object rval, 
                                bits16 child_count, object reader ):
        cdef bits32 start_chrom_ix, start_base, end_chrom_ix, end_base
        cdef bits64 offset
        cdef bits64 size
        for i from 0 <= i < child_count:
            start_chrom_ix = reader.read_uint32()
            start_base = reader.read_uint32()
            end_chrom_ix = reader.read_uint32()
            end_base = reader.read_uint32()
            offset = reader.read_uint64()
            size = reader.read_uint64()
            if overlaps( chrom_ix, start, end, start_chrom_ix, start_base, end_chrom_ix, end_base ):
                rval.append( ( offset, size ) )

    def r_find_overlapping_parent( self, int level, bits32 chrom_ix, bits32 start, bits32 end, object rval, 
                                  bits16 child_count, object reader ):
        # Read and cache offsets for all children to avoid excessive seeking
        ## cdef bits32 start_chrom_ix[child_count], start_base[child_count], end_chrom_ix[child_count], end_base[child_count]
        ## cdef bits64 offset[child_count]
        start_chrom_ix = []; start_base = []; end_chrom_ix = []; end_base = []
        offset = []
        for i from 0 <= i < child_count:
            ## start_chrom_ix[i] = reader.read_bits32()
            ## start_base[i] = reader.read_bits32()
            ## end_chrom_ix[i] = reader.read_bits32()
            ## end_base[i] = reader.read_bits32()
            ## offset[i] = reader.read_bits64()
            start_chrom_ix.append( reader.read_uint32() )
            start_base.append( reader.read_uint32() )
            end_chrom_ix.append( reader.read_uint32() )
            end_base.append( reader.read_uint32() )
            offset.append( reader.read_uint64() )
        # Now recurse
        for i from 0 <= i < child_count:
            if overlaps( chrom_ix, start, end, start_chrom_ix[i], start_base[i], end_chrom_ix[i], end_base[i] ):
                self.r_find_overlapping( level + 1, offset[i], chrom_ix, start, end, rval, reader )

    def find_overlapping_blocks( self, bits32 chrom_ix, bits32 start, bits32 end ):
        rval = []
        self.r_find_overlapping( 0, self.root_offset, chrom_ix, start, end, rval, self.reader )
        return rval
