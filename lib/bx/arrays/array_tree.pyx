from __future__ import division

__all__ = [ 'ArrayTree', 'FileArrayTreeDict', 'array_tree_dict_from_reader' ]

import numpy
from numpy import *
cimport numpy

cimport bx.arrays.wiggle

from bx.misc.binary_file import BinaryFileWriter, BinaryFileReader
from bx.misc.cdb import FileCDBDict

"""
Classes for storing binary data on disk in a tree structure that allows for
efficient sparse storage (when the data occurs in contiguous blocks), fast
access to a specific block of data, and fast access to summaries at different
resolutions.

On disk format
--------------

Blocks are stored contiguously on disk in level-order. Contents should always be
network byte order (big endian), however this implementation will byte-swap when
reading  if necessary. File contents:

- magic:      uint32
- version:    unit32
- array size: uint32
- block size: uint32
- array type: 4 chars (numpy typecode, currently only simple types represented by one char are supported)

- Internal nodes in level order
    - Summary
        - count of valid values in each subtree : sizeof( dtype ) * block_size
        - frequencies: sizeof ( int32 ) * block_size
        - min of valid values in each subtree : sizeof( dtype ) * block_size
        - max of valid values in each subtree : sizeof( dtype ) * block_size
        - sum of valid values in each subtree : sizeof( dtype ) * block_size
        - sum of squares of valid values in each subtree : sizeof( dtype ) * block_size
    - File offsets of each child node: uint64 * block_size
  
- Leaf nodes
    - data points: sizeof( dtype ) * block_size
    
- Version 1 reads version 0 and version 1

"""

## Enhancement ideas:
## 
##   - Write markers of the number of blocks skipped between blocks. This would
##     allow fast striding across data or summaries (use the indexes to get to
##     the start of a block, then just read straight through). Would this help?
## 
##   - Compression for blocks?

MAGIC = 0x310ec7dc
VERSION = 1
NUM_SUMMARY_ARRAYS = 6

def array_tree_dict_from_reader( reader, sizes, default_size=2147483647, block_size=1000, no_leaves=False ):
    # Create empty array trees
    rval = {}
    ## for key, size in sizes.iteritems():
    ##    rval[ key ] = ArrayTree( size, 1000 )
    # Fill
    last_chrom = None
    last_array_tree = None
    for chrom, start, end, _, val in reader:
        if chrom != last_chrom:
            if chrom not in rval:
                rval[chrom] = ArrayTree( sizes.get( chrom, default_size ), block_size, no_leaves=no_leaves )
            last_array_tree = rval[chrom]
        last_array_tree.set_range( start, end, val )
    return rval
                

cdef class FileArrayTreeDict:
    """
    Access to a file containing multiple array trees indexed by a string key.
    """
    cdef object io
    cdef object cdb_dict
    def __init__( self, file ):
        self.io = io = BinaryFileReader( file, MAGIC )
        assert (0 <= io.read_uint32() <= 1) # Check for version 0 or 1
        self.cdb_dict = FileCDBDict( file, is_little_endian=io.is_little_endian )
    def __getitem__( self, key ):
        offset = self.cdb_dict[key]
        offset = self.io.unpack( "L", offset.encode() )[0]
        self.io.seek( offset )
        return FileArrayTree( self.io.file, self.io.is_little_endian )
    
    @classmethod
    def dict_to_file( Class, dict, file, is_little_endian=True, no_leaves=False ):
        """
        Writes a dictionary of array trees to a file that can then be
        read efficiently using this class.
        """
        io = BinaryFileWriter( file, is_little_endian=is_little_endian )
        # Write magic number and version
        io.write_uint32( MAGIC )
        io.write_uint32( VERSION )
        # Write cdb index with fake values just to fill space
        cdb_dict = {}
        for key in dict.iterkeys():
            cdb_dict[ key ] = io.pack( "L", 0 )
        cdb_offset = io.tell()
        FileCDBDict.to_file( cdb_dict, file, is_little_endian=is_little_endian )
        # Write each tree and save offset
        for key, value in dict.iteritems():
            offset = io.tell()
            cdb_dict[ key ] = io.pack( "L", offset )
            value.to_file( file, is_little_endian=is_little_endian, no_leaves=no_leaves )
        # Go back and write the index again
        io.seek( cdb_offset )
        FileCDBDict.to_file( cdb_dict, file, is_little_endian=is_little_endian )
        
cdef class FileArrayTree:
    """
    Wrapper for ArrayTree stored in file that reads as little as possible
    """
    cdef public int max
    cdef public int block_size
    cdef public object dtype
    cdef public int levels
    cdef public int offset
    cdef public int root_offset
    cdef object io
    
    def __init__( self, file, is_little_endian=True ):
        self.io = BinaryFileReader( file, is_little_endian=is_little_endian )
        self.offset = self.io.tell()
        # Read basic info about the tree
        self.max = self.io.read_uint32()
        self.block_size = self.io.read_uint32()
        # Read dtype and canonicalize
        dt = self.io.read( 1 )
        self.dtype = numpy.dtype( dt )
        self.io.skip( 3 )
        # How many levels are needed to cover the entire range?
        self.levels = 0
        while ( <long long> self.block_size ) ** ( self.levels + 1 ) < self.max:
            self.levels += 1
        # Not yet dealing with the case where the root is a Leaf
        assert self.levels > 0, "max < block_size not yet handled"
        # Save offset of root
        self.root_offset = self.io.tell()
        
    def __getitem__( self, index ):
        min = self.r_seek_to_node( index, 0, self.root_offset, self.levels, 0 )
        if min < 0:
            return nan
        self.io.skip( self.dtype.itemsize * ( index - min ) )
        return self.io.read_raw_array( self.dtype, 1 )[0]
        
    def get_summary( self, index, level ):
        if level <= 0 or level > self.levels:
            raise ValueError, "level must be <= self.levels"
        if self.r_seek_to_node( index, 0, self.root_offset, self.levels, level ) < 0:
            return None
        # Read summary arrays
        s = Summary()
        s.counts = self.io.read_raw_array( self.dtype, self.block_size )
        s.frequencies = self.io.read_raw_array( self.dtype, self.block_size )
        s.sums = self.io.read_raw_array( self.dtype, self.block_size )
        s.mins = self.io.read_raw_array( self.dtype, self.block_size)
        s.maxs = self.io.read_raw_array( self.dtype, self.block_size )
        s.sumsquares = self.io.read_raw_array( self.dtype, self.block_size )
        return s
        
    def get_leaf( self, index ):
        if self.r_seek_to_node( index, 0, self.root_offset, self.levels, 0 ) < 0:
            return []
        return self.io.read_raw_array( self.dtype, self.block_size )
        
    cdef int r_seek_to_node( self, int index, int min, long long offset, int level, int desired_level ):
        """
        Seek to the start of the node at `desired_level` that contains `index`.
        Returns the minimum value represented in that node.
        """
        cdef int child_size, bin_index, child_min
        self.io.seek( offset )
        if level > desired_level:
            child_size = self.block_size ** level
            bin_index = ( index - min ) // ( child_size ) 
            child_min = min + ( bin_index * child_size )
            # Skip summary arrays -- # arrays * itemsize * block_size
            self.io.skip( NUM_SUMMARY_ARRAYS * self.dtype.itemsize * self.block_size )
            # Skip to offset of correct child -- offsets are 8 bytes
            self.io.skip( 8 * bin_index )
            # Read offset of child
            child_offset = self.io.read_uint64()
            # print "co: %s\tbi: %s\tcm: %s\n" % (child_offset, bin_index, child_min)
            if child_offset == 0:
                return -1
            return self.r_seek_to_node( index, child_min, child_offset, level - 1, desired_level )
        else:
            # The file pointer is at the start of the desired node, do nothing
            return min

cdef class Summary:
    """
    Summary for a non-leaf level of the tree, contains arrays of the min, max,
    valid count, sum, and sum-of-squares for each child.
    """
    cdef public object counts
    cdef public object frequencies
    cdef public object mins
    cdef public object maxs
    cdef public object sums
    cdef public object sumsquares

cdef class ArrayTreeNode
cdef class ArrayTreeLeaf

cdef class ArrayTree:
    """
    Stores a sparse array of data as a tree.
    
    An array of `self.max` values is stored in a tree in which each leaf
    contains `self.block_size` values and each internal node contains
    `self.block_size` children.
    
    Entirely empty subtrees are not stored. Thus, the storage is efficient for
    data that is block sparse -- having contiguous chunks of `self.block_size` or
    larger data. Currently it is not efficient if the data is strided (e.g.
    one or two data points in every interval of length `self.block_size`).
    
    Internal nodes store `Summary` instances for their subtrees. 
    """

    cdef public int max
    cdef public int block_size
    cdef public object dtype
    cdef public int levels
    cdef public int no_leaves
    cdef public ArrayTreeNode root

    def __init__( self, int max, int block_size, dtype=float32, no_leaves=False ):
        """
        Create a new array tree of size `max` 
        """
        self.max = max
        self.block_size = block_size
        self.no_leaves = no_leaves
        # Force the dtype argument to its canonical dtype object
        self.dtype = numpy.dtype( dtype )
        # How many levels are needed to cover the entire range?
        self.levels = 0
        while ( <long long> self.block_size ) ** ( self.levels + 1 ) < self.max:
            self.levels += 1
        # Not yet dealing with the case where the root is a Leaf
        assert self.levels > 0, "max < block_size not yet handled"
        # Create the root node`
        self.root = ArrayTreeNode( self, 0, max, block_size, self.levels )
        
    def __setitem__( self, int index, value ):
        self.root.set( index, value )
        
    def set_range( self, int start, int end, value ):
        for i from start <= i < end:
            self.root.set( i, value )
        
    def __getitem__( self, int index ):
        return self.root.get( index )
        
    def to_file( self, f, is_little_endian=True, no_leaves=False ):
        io = BinaryFileWriter( f, is_little_endian=is_little_endian )
        ## io.write_uint32( VERSION )
        io.write_uint32( self.max )
        io.write_uint32( self.block_size )
        io.write( self.dtype.char )
        io.write( "\0\0\0" )
        # Data pass, level order
        if no_leaves:
            bottom_level = 0
        else:
            bottom_level = -1
        for level in range( self.levels, bottom_level, -1 ):
            self.root.to_file_data_pass( io, level )
        # Offset pass to fix up indexes
        self.root.to_file_offset_pass( io )
        
    @classmethod
    def from_file( Class, f, is_little_endian=True ):        
        io = BinaryFileReader( f, is_little_endian=is_little_endian )
        ## assert io.read_uint32() == VERSION
        max = io.read_uint32()
        block_size = io.read_uint32()
        dt = io.read( 1 )
        io.read( 3 )
        tree = Class( max, block_size, dt )
        tree.root.from_file( io )
        return tree
    
    @classmethod
    def from_sequence( Class, s, block_size=1000 ):
        """
        Build an ArrayTree from a sequence like object (must have at least
        length and getitem).
        """
        tree = Class( len( s ), block_size )
        for i in range( len( s ) ):
            tree[i] = s[i]
        return tree
        
cdef class ArrayTreeNode:
    """
    Internal node of an ArrayTree. Contains summary data and pointers to
    subtrees.
    """

    cdef ArrayTree tree
    cdef int min
    cdef int max
    cdef int block_size
    cdef int level
    cdef int child_size
    cdef object children
    cdef public Summary summary
    cdef public long start_offset

    def __init__( self, ArrayTree tree, int min, int max, int block_size, int level ):
        self.tree = tree
        self.min = min
        self.max = max
        self.block_size = block_size
        self.level = level
        # Each of my children represents block_size ** level values
        self.child_size = self.block_size ** self.level        
        self.children = [None] * self.block_size
        self.summary = None
        self.start_offset = 0
        
    cdef inline init_bin( self, int index ):
        cdef int min = self.min + ( index * self.child_size )
        cdef int max = min + self.child_size
        if self.level == 1:
            self.children[ index ] = ArrayTreeLeaf( self.tree, min, max )
        else:
            self.children[ index ] = ArrayTreeNode( self.tree, min, max, self.block_size, self.level - 1 )
            
    def set( self, int index, value ):
        cdef int bin_index = ( index - self.min ) // ( self.child_size )
        if self.children[ bin_index ] is None:
            self.init_bin( bin_index )
        self.children[ bin_index ].set( index, value )
        
    def get( self, int index ):
        cdef int bin_index = ( index - self.min ) // ( self.child_size ) 
        if self.children[ bin_index ] is None:
            return nan
        else:
            return self.children[ bin_index ].get( index )
            
    cpdef build_summary( self ):
        """
        Build summary of children. 
        """
        counts = empty( self.tree.block_size, self.tree.dtype )
        frequencies = empty( self.tree.block_size, self.tree.dtype )
        mins = empty( self.tree.block_size, self.tree.dtype )
        maxs = empty( self.tree.block_size, self.tree.dtype )
        sums = empty( self.tree.block_size, self.tree.dtype )
        sumsquares = empty( self.tree.block_size, self.tree.dtype )
        for i in range( len( self.children ) ):
            if self.children[i]:
                if self.level == 1:
                    v = self.children[i].values
                    counts[i] = sum( ~isnan( v ) )
                    frequencies[i] = self.children[i].frequency
                    mins[i] = nanmin( v )
                    maxs[i] = nanmax( v )
                    sums[i] = nansum( v )
                    sumsquares[i] = nansum( v ** 2 )
                else:
                    c = self.children[i]
                    c.build_summary()
                    counts[i] = sum( c.summary.counts )
                    frequencies[i] = sum( c.summary.frequencies )
                    mins[i] = nanmin( c.summary.mins )
                    maxs[i] = nanmax( c.summary.maxs )
                    sums[i] = nansum( c.summary.sums )
                    sumsquares[i] = nansum( c.summary.sumsquares )
            else:
                counts[i] = 0
                frequencies[i] = 0
                mins[i] = nan
                maxs[i] = nan
                sums[i] = nan
                sumsquares[i] = nan
        s = Summary()
        s.counts = counts
        s.frequencies = frequencies
        s.mins = mins
        s.maxs = maxs
        s.sums = sums
        s.sumsquares = sumsquares
        self.summary = s
        
    def to_file_data_pass( self, io, level ):
        """
        First pass of writing to file, writes data and saves position of block.
        """
        assert self.summary, "Writing without summaries is currently not supported"
        # If we are at the current level being written, write a block
        if self.level == level:
            # Save file offset where this block starts
            self.start_offset = io.tell()
            # Write out summary data
            io.write_raw_array( self.summary.counts )
            io.write_raw_array( self.summary.frequencies )
            io.write_raw_array( self.summary.sums )
            io.write_raw_array( self.summary.mins )
            io.write_raw_array( self.summary.maxs )
            io.write_raw_array( self.summary.sumsquares )
            # Skip enough room for child offsets (block_size children * 64bits)
            io.skip( self.tree.block_size * 8 )
        # Must be writing a lower level, so recurse
        else:
            # Write all non-empty children
            for i in range( len( self.children ) ):
                if self.children[i] is not None:
                    self.children[i].to_file_data_pass( io, level )
                
    def to_file_offset_pass( self, io ):
        """
        Second pass of writing to file, seek to appropriate position and write
        offsets of children.
        """
        # Seek to location of child offfsets (skip over # summary arrays)
        skip_amount = NUM_SUMMARY_ARRAYS * self.tree.dtype.itemsize * self.block_size
        io.seek( self.start_offset + skip_amount )
        # Write the file offset of each child into the index
        for child in self.children:
            if child is None:
                io.write_uint64( 0 )
            else:
                io.write_uint64( child.start_offset )
        # Recursively write offsets in child nodes
        for child in self.children:
            if child is not None:
                child.to_file_offset_pass( io )
    
    def from_file( self, io ):
        """
        Load entire summary and all children into memory.
        """
        dtype = self.tree.dtype
        block_size = self.tree.block_size
        # Read summary arrays
        s = Summary()
        s.counts = io.read_raw_array( dtype, block_size )
        s.frequencies = io.read_raw_array( int32, block_size )
        s.sums = io.read_raw_array( dtype, block_size )
        s.mins = io.read_raw_array( dtype, block_size)
        s.maxs = io.read_raw_array( dtype, block_size )
        s.sumsquares = io.read_raw_array( dtype, block_size )
        self.summary = s
        # Read offset of all children
        child_offsets = [ io.read_uint64() for i in range( block_size ) ]
        for i in range( block_size ):
            if child_offsets[i] > 0:
                self.init_bin( i )
                io.seek( child_offsets[i] )
                self.children[i].from_file( io )
                
    def get_from_file( self, io, index ):
        cdef int bin_index = ( index - self.min ) //( self.child_size ) 
        if self.children[ bin_index ] is None:
            return nan
        else:
            return self.children[ bin_index ].get( index )
        
cdef class ArrayTreeLeaf:
    """
    Leaf node of an ArrayTree, contains data values.
    """
    
    cdef ArrayTree tree
    cdef int min
    cdef int max
    cdef public int frequency
    cdef public numpy.ndarray values
    cdef public long start_offset
    
    def __init__( self, ArrayTree tree, int min, int max ):
        self.tree = tree
        self.min = min
        self.max = max
        self.frequency = 0
        self.values = empty( max - min, self.tree.dtype )
        self.values[:] = nan
        self.start_offset = 0
        
    def set( self, index, value ):
        self.frequency += 1
        self.values[ index - self.min ] = value
        
    def get( self, index ):
        return self.values[ index - self.min ]
        
    def to_file_data_pass( self, io, level ):
        assert level == 0
        self.start_offset = io.tell()
        io.write_raw_array( self.values )
    
    def to_file_offset_pass( self, io ):
        pass
        
    def from_file( self, io ):
        self.values = io.read_raw_array( self.tree.dtype, self.tree.block_size )
