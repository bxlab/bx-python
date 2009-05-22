from __future__ import division

__all__ = [ 'ArrayTree' ]

import numpy
from numpy import *
cimport numpy

from bx.misc.binary_file import BinaryFileWriter, BinaryFileReader

"""
Classes for storing binary data on disk in a tree structure that allows for
efficient sparse storage (when the data occurs in contiguous blocks), fast
access to a specific block of data, and fast access to summaries at different
resolutions.

NOTE: This is still a work in progress, the file format for the first version
      is not yet finalized, once it is the version will be incremented to 1.

On disk format
--------------

Blocks are stored contigously on disk in level-order. Contents should always be
network byte order (big endian), however this implementation will byte-swap when
reading  if neccesary. File contents:

- magic:      uint32
- version:    unit32
- array size: uint32
- block size: uint32
- array type: 4 chars (numpy typecode, currently only simple types represented by one char are supported)

- Internal nodes in level order
  - Summary
      - count of valid values in each subtree : sizeof( dtype ) * block_size
      - min of valid values in each subtree : sizeof( dtype ) * block_size
      - max of valid values in each subtree : sizeof( dtype ) * block_size
      - sum of valid values in each subtree : sizeof( dtype ) * block_size
      - sum of squares of  of valid values in each subtree : sizeof( dtype ) * block_size
  - File offsets of each child node: uint64 * block_size
  
- Leaf nodes
  - data points: sizeof( dtype ) * block_size

"""

## Enhancement ideas:
## 
##   - Write markers of the number of blocks skipped between blocks. This would
##     allow fast striding across data or summaries (use the indexes to get to
##     the start of a block, then just read straight through). Would this help?
## 
##   - Compression for blocks?

MAGIC = 0x310ec7dc
VERSION = 0

cdef class Summary( object ):
    """
    Summary for a non-leaf level of the tree, contains arrays of the min, max,
    valid count, sum, and sum-of-squares for each child.
    """
    cdef public object counts
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

    cdef int max
    cdef int block_size
    cdef object dtype
    cdef int levels
    cdef public ArrayTreeNode root

    def __init__( self, int max, int block_size, dtype=float32 ):
        """
        Create a new array tree of size `max` 
        """
        self.max = max
        self.block_size = block_size
        # Force the dtype argument to its canonical dtype object
        self.dtype = numpy.dtype( dtype )
        # How many levels are needed to cover the entire range?
        self.levels = 0
        while self.block_size ** ( self.levels + 1 ) < self.max:
            self.levels += 1
        # Not yet dealing with the case where the root is a Leaf
        assert self.levels > 0, "max < block_size not yet handled"
        # Create the root node
        self.root = ArrayTreeNode( self, 0, max, block_size, self.levels )
        
    def __setitem__( self, int index, value ):
        self.root.set( index, value )
        
    def __getitem__( self, int index ):
        return self.root.get( index )
        
    def to_file( self, f ):
        io = BinaryFileWriter( f, magic=MAGIC )
        io.write_uint32( VERSION )
        io.write_uint32( self.max )
        io.write_uint32( self.block_size )
        io.write( self.dtype.char )
        io.write( "\0\0\0" )
        # Data pass, level order
        for level in range( self.levels, -1, -1 ):
            self.root.to_file_data_pass( io, level )
        # Offset pass to fix up indexes
        self.root.to_file_offset_pass( io )
        
    @classmethod
    def from_file( Class, f ):        
        io = BinaryFileReader( f, magic=MAGIC )
        assert io.read_uint32() == VERSION
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
        cdef int bin_index = ( index - self.min ) //( self.child_size ) 
        if self.children[ bin_index ] is None:
            return nan
        else:
            return self.children[ bin_index ].get( index )
            
    cpdef build_summary( self ):
        """
        Build summary of children. 
        """
        counts = empty( self.tree.block_size, self.tree.dtype )
        mins = empty( self.tree.block_size, self.tree.dtype )
        maxs = empty( self.tree.block_size, self.tree.dtype )
        sums = empty( self.tree.block_size, self.tree.dtype )
        sumsquares = empty( self.tree.block_size, self.tree.dtype )
        for i in range( len( self.children ) ):
            if self.children[i]:
                if self.level == 1:
                    v = self.children[i].values
                    counts[i] = sum( ~isnan( v ) )
                    mins[i] = nanmin( v )
                    maxs[i] = nanmax( v )
                    sums[i] = nansum( v )
                    sumsquares[i] = nansum( v ** 2 )
                else:
                    c = self.children[i]
                    c.build_summary()
                    counts[i] = sum( c.summary.counts )
                    mins[i] = nanmin( c.summary.mins )
                    maxs[i] = nanmax( c.summary.maxs )
                    sums[i] = nansum( c.summary.sums )
                    sumsquares[i] = nansum( c.summary.sumsquares )
            else:
                counts[i] = 0
                mins[i] = nan
                maxs[i] = nan
                sums[i] = nan
                sumsquares[i] = nan
        s = Summary()
        s.counts = counts
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
        # Seek to location of child offfsets (skip over 5 summary arrays)
        skip_amount = 5 * self.tree.dtype.itemsize * self.block_size
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
        
cdef class ArrayTreeLeaf:
    """
    Leaf node of an ArrayTree, contains data values.
    """
    
    cdef ArrayTree tree
    cdef int min
    cdef int max
    cdef public numpy.ndarray values
    cdef public long start_offset
    
    def __init__( self, ArrayTree tree, int min, int max ):
        self.tree = tree
        self.min = min
        self.max = max
        self.values = empty( max - min, self.tree.dtype )
        self.values[:] = nan
        self.start_offset = 0
        
    def set( self, index, value ):
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
