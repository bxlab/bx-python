from align import *

import itertools
import interval_index_file

# Tools for dealing with pairwise alignments in LAV format

class MultiIndexed( object ):
    """Similar to 'indexed' but wraps more than one lav_file"""
    def __init__( self, lav_filenames, keep_open=False ):
        self.indexes = [ Indexed( lav_file, lav_file + ".index" ) for lav_file in lav_filenames ]
    def get( self, src, start, end ):
        blocks = []
        for index in self.indexes: blocks += index.get( src, start, end )
        return blocks


class Indexed( object ):
    """Indexed access to a lav using overlap queries, requires an index file"""

    def __init__( self, lav_filename, index_filename=None, keep_open=False, species_to_lengths=None ):
        if index_filename is None: index_filename = lav_filename + ".index"
        self.indexes = interval_index_file.Indexes( filename=index_filename )
        self.lav_filename = lav_filename
        self.species_to_lengths = species_to_lengths
        if keep_open:
            self.f = open( lav_filename )
        else:
            self.f = None

    def get( self, src, start, end ):
        intersections = self.indexes.find( src, start, end )
        return itertools.imap( self.get_lav_at_offset, [ val for start, end, val in intersections ] )

    def get_lav_at_offset( self, offset ):
        if self.f:
            self.f.seek( offset )
            return read_next_lav( self.f, self.species_to_lengths )
        else:
            f = open( self.lav_filename )
            try:
                f.seek( offset )
                return read_next_lav( f, self.species_to_lengths )
            finally:
                f.close()

class Reader( object ):
    """Iterate over all lav blocks in a file in order"""
    
    def __init__( self, file, species_to_lengths=None ):
        self.file = file
        self.species_to_lengths = species_to_lengths
        self.attributes = {}

    def next( self ):
        return read_next_lav( self.file, self.species_to_lengths )

    def __iter__( self ):
        return ReaderIter( self )

    def close( self ):
        self.file.close()

class ReaderIter( object ):
    def __init__( self, reader ):
        self.reader = reader
    def __iter__( self ): 
        return self
    def next( self ):
        v = self.reader.next()
        if not v: raise StopIteration
        return v

class Writer( object ):

    def __init__( self, file, attributes={} ):
        self.file = file

    def write( self, alignment ):
        if (len(alignment.components) != 2):
            raise "%d-component alignment is not compatible with lav" % \
                   len(alignment.components)

        raise "lav writing is not yet supported"

    def close( self ):
        self.file.close()

# ---- Helper methods ---------------------------------------------------------

def read_next_lav( file, species_to_lengths=None ):
    raise "lav reading is not yet supported"

def readline( file, skip_blank=False ):
    """Read a line from provided file, skipping any blank or comment lines"""
    while 1:
        line = file.readline()
        if not line: return None
        if line[0] != '#' and not ( skip_blank and line.isspace() ):
            return line

