from __future__ import division

import sys
from bx_extras.lrucache import LRUCache
from cStringIO import StringIO

DEFAULT_CACHE_SIZE=10
DEFAULT_BLOCK_SIZE=1024*1024*2

class FileCache( object ):
    """
    Wrapper for a file that cache blocks of data in memory. 
    
    **NOTE:** this is currently an incomplete file-like object, it only
    supports seek, tell, and readline (plus iteration). Reading bytes is
    currently not implemented.
    """
    def __init__( self, file, size, cache_size=DEFAULT_CACHE_SIZE, 
                                    block_size=DEFAULT_BLOCK_SIZE ):
        """
        Create a new `FileCache` wrapping the file-like object `file` that
        has total size `size` and caching blocks of size `block_size`.
        """
        self.file = file
        self.size = size
        self.cache_size = cache_size
        self.block_size = block_size
        # Setup the cache
        self.nblocks = ( self.size // self.block_size ) + 1
        self.cache = LRUCache( self.cache_size )
        # Position in file
        self.dirty = True
        self.at_eof = False
        self.file_pos = 0
        self.current_block_index = -1
        self.current_block = None
    def fix_dirty( self ):
        chunk, offset = self.get_block_and_offset( self.file_pos )
        if self.current_block_index != chunk:
            self.current_block = StringIO( self.load_block( chunk ) )
            self.current_block.read( offset )
            self.current_block_index = chunk
        else:
            self.current_block.seek( offset )
        self.dirty = False
    def get_block_and_offset( self, index ):
        return int( index // self.block_size ), int( index % self.block_size )
    def load_block( self, index ):
        if index in self.cache:
            return self.cache[index]
        else:
            real_offset = index * self.block_size
            self.file.seek( real_offset )
            block = self.file.read( self.block_size )
            self.cache[index] = block
            return block
    def seek( self, offset, whence=0 ):
        """
        Move the file pointer to a particular offset.
        """
        # Determine absolute target position
        if whence == 0:
            target_pos = offset
        elif whence == 1:
            target_pos = self.file_pos + offset
        elif whence == 2:
            target_pos = self.size - offset
        else:
            raise Exception( "Invalid `whence` argument: %r", whence )
        # Check if this is a noop
        if target_pos == self.file_pos:
            return    
        # Verify it is valid
        assert 0 <= target_pos < self.size, "Attempt to seek outside file"
        # Move the position
        self.file_pos = target_pos
        # Mark as dirty, the next time a read is done we need to actually
        # move the position in the bzip2 file
        self.dirty = True
    def readline( self ):
        if self.dirty:
            self.fix_dirty()
        if self.at_eof:
            return ""
        rval = []
        while 1:
            line = self.current_block.readline()
            rval.append( line )
            if len( line ) > 0 and line[-1] == '\n':
                break
            elif self.current_block_index == self.nblocks - 1:
                self.at_eof = True
                break
            else:
                self.current_block_index += 1
                self.current_block = StringIO( self.load_block( self.current_block_index ) )      
        return "".join( rval )     
    def next( self ):
        line = self.readline()
        if line == "":
            raise StopIteration
    def __iter__( self ):
        return self
    def close( self ):
        self.file.close()