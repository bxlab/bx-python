from bx.align import *

import itertools
from bx import interval_index_file

# Tools for dealing with multiple alignments in MAF format

class MultiIndexed( object ):
    """Similar to 'indexed' but wraps more than one maf_file"""
    def __init__( self, maf_filenames, keep_open=False ):
        self.indexes = [ Indexed( maf_file, maf_file + ".index" ) for maf_file in maf_filenames ]
    def get( self, src, start, end ):
        blocks = []
        for index in self.indexes: blocks.extend( index.get( src, start, end ) )
        return blocks
    
class Indexed( object ):
    """Indexed access to a maf using overlap queries, requires an index file"""

    def __init__( self, maf_filename, index_filename=None, keep_open=False, species_to_lengths=None ):
        if index_filename is None: index_filename = maf_filename + ".index"
        self.indexes = interval_index_file.Indexes( filename=index_filename )
        self.maf_filename = maf_filename
        self.species_to_lengths = species_to_lengths
        if keep_open: 
            self.f = open( maf_filename )
        else:
            self.f = None

    def get( self, src, start, end ):
        intersections = self.indexes.find( src, start, end )
        return map( self.get_maf_at_offset, [ val for start, end, val in intersections ] )

    def get_maf_at_offset( self, offset ):
        if self.f:
            self.f.seek( offset )
            return read_next_maf( self.f, self.species_to_lengths )
        else:
            f = open( self.maf_filename )
            try:
                f.seek( offset )
                return read_next_maf( f, self.species_to_lengths ) 
            finally:
                f.close()
            
class Reader( object ):
    """Iterate over all maf blocks in a file in order"""
    
    def __init__( self, file, species_to_lengths=None ):
        self.file = file
        self.species_to_lengths = species_to_lengths
        # Read and verify maf header, store any attributes
        fields = self.file.readline().split()
        if fields[0] != '##maf': raise "File does not have MAF header"
        self.attributes = parse_attributes( fields[1:] )

    def next( self ):
        return read_next_maf( self.file, self.species_to_lengths )

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
        # Write header, Webb's maf code wants version first, we accomodate
        if not attributes.has_key('version'): attributes['version'] = 1
        self.file.write( "##maf version=%s" % attributes['version'] )
        for key in attributes: 
            if key == 'version': continue
            self.file.writelines( " %s=%s" % ( key, attributes[key] ) )
        self.file.write( "\n" )

    def write( self, alignment ):
        self.file.write( "a score=" + str( alignment.score ) )
        for key in alignment.attributes:
            self.file.write( " %s=%s" % ( key, alignment.attributes[key] ) )
        self.file.write( "\n" )
        # Components
        rows = []
        for c in alignment.components:
            rows.append( ( "s", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.text ) )
        self.file.write( format_tabular( rows, "llrrrrr" ) )
        self.file.write( "\n" )

    def close( self ):
        self.file.close()

# ---- Helper methods ---------------------------------------------------------

def read_next_maf( file, species_to_lengths=None ):
    alignment = Alignment(species_to_lengths=species_to_lengths)
    # Attributes line
    line = readline( file, skip_blank=True )
    if not line: return None
    fields = line.split() 
    if fields[0] != 'a': raise "Expected 'a ...' line"
    alignment.attributes = parse_attributes( fields[1:] )
    alignment.score = alignment.attributes['score']
    del alignment.attributes['score']
    # Sequence lines
    while 1:
        line = readline( file )
        # EOF or Blank line terminates alignment components
        if not line or line.isspace(): break
        if line.isspace(): break 
        # Verify
        fields = line.split()
        if fields[0] != 's': raise "Expected 's ...' line"
        # Parse 
        component = Component()
        component.src = fields[1]
        component.start = int( fields[2] )
        component.size = int( fields[3] )
        component.strand = fields[4]
        component.src_size = int( fields[5] )
        if len(fields) > 6: component.text = fields[6].strip()
        # Add to set
        alignment.add_component( component )
    return alignment

def readline( file, skip_blank=False ):
    """Read a line from provided file, skipping any blank or comment lines"""
    while 1:
        line = file.readline()
        if not line: return None 
        if line[0] != '#' and not ( skip_blank and line.isspace() ):
            return line

def parse_attributes( fields ):
    """Parse list of key=value strings into a dict"""
    attributes = {}
    for field in fields:
        pair = field.split( '=' )
        attributes[ pair[0] ] = pair[1]
    return attributes

def format_tabular( rows, align=None ):
    if len( rows ) == 0: return ""
    lengths = [ len( col ) for col in rows[ 0 ] ]
    for row in rows[1:]:
        for i in range( 0, len( row ) ):
            lengths[ i ] = max( lengths[ i ], len( row[ i ] ) )
    rval = ""
    for row in rows:
        for i in range( 0, len( row ) ):
            if align and align[ i ] == "l":
                rval += row[ i ].ljust( lengths[ i ] )
            else:
                rval += row[ i ].rjust( lengths[ i ] )
            rval += " "
        rval += "\n"
    return rval
        
