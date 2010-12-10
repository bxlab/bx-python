from bx.misc.binary_file import BinaryFileReader

DEF bpt_sig = 0x78CA8C91

# bptFileHeaderSize = 32
# bptBlockHeaderSize = 4

cdef class BPTFile:
    """
    On disk B+ tree compatible with Jim Kent's bPlusTree.c
    """

    def __init__( self, file=None ):
        if file is not None:
            self.attach( file )

    def attach( self, file ):
        """
        Attach to an open file
        """
        self.file = file
        self.reader = reader = BinaryFileReader( file, bpt_sig )
        self.is_byteswapped = self.reader.byteswap_needed
        # Read header stuff
        self.block_size = reader.read_uint32()
        self.key_size = reader.read_uint32()
        self.value_size = reader.read_uint32()
        self.item_count = reader.read_uint64()
        reader.skip( 8 )
        self.root_offset = reader.tell()

    def r_find( self, bits64 block_start, key ):
        """
        Recursively seek the value matching key under the subtree starting
        at file offset `block_start`
        """
        cdef UBYTE is_leaf
        cdef bits16 child_count
        cdef bits64 offset
        self.reader.seek( block_start )
        # Block header
        is_leaf = self.reader.read_uint8()
        self.reader.read_uint8()
        child_count = self.reader.read_uint16()
        if is_leaf:
            for i from 0 <= i < child_count:
                node_key = self.reader.read( self.key_size )
                node_value = self.reader.read( self.value_size )
                if node_key == key:
                    return node_value
            return None
        else:
            # Read and discard first key, store offset
            self.reader.read( self.key_size )
            offset = self.reader.read_bits64()
            # Loop until correct subtree is found
            for i from 0 <= i < child_count:
                node_key = self.reader.read( self.key_size )
                if node_key > key:
                    break
                offset = self.reader.read_bits64()
            return self.r_find( offset, key )

    def find( self, key ):
        """
        Find the value matching `key` (a bytestring). Returns the matching
        value as a bytestring if found, or None
        """
        # Key is greater than key_size, must not be a match
        if len(key) > self.key_size:
            return None
        # Key is less than key_size, right pad with 0 bytes
        if len(key) < self.key_size:
            key += ( '\0' * ( self.key_size - len(key) ) )
        # Call the recursive finder
        return self.r_find( self.root_offset, key )
