"""
Support for the `MAF`_ multiple sequence alignment format used by `multiz`_.

.. _MAF: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
.. _multiz: http://www.bx.psu.edu/miller_lab/
"""

from bx.align import *

from StringIO import StringIO
import os

import itertools
from bx import interval_index_file

from bx.misc.seekbzip2 import SeekableBzip2File

MAF_INVERSE_STATUS = 'V'
MAF_INSERT_STATUS = 'I'
MAF_CONTIG_STATUS = 'C'
MAF_CONTIG_NESTED_STATUS = 'c'
MAF_NEW_STATUS = 'N'
MAF_NEW_NESTED_STATUS = 'n'
MAF_MAYBE_NEW_STATUS = 'S'
MAF_MAYBE_NEW_NESTED_STATUS = 's'
MAF_MISSING_STATUS = 'M'

class MAFIndexedAccess( interval_index_file.AbstractIndexedAccess ):
    """
    Indexed access to a MAF file.
    """
    def read_at_current_offset( self, file, **kwargs ):
        """
        Read the MAF block at the current position in `file` and return an
        instance of `Alignment`.
        """
        return read_next_maf( file, **kwargs )

class MAFMultiIndexedAccess( interval_index_file.AbstractMultiIndexedAccess ):
    """
    Indexed access to multiple MAF files.
    """
    indexed_access_class = MAFIndexedAccess
      
Indexed = MAFIndexedAccess
"""Deprecated: `MAFIndexedAccess` is also available under the name `Indexed`."""

MultiIndexed = MAFMultiIndexedAccess
"""Deprecated: `MAFMultiIndexedAccess` is also available under the name `MultiIndexed`."""

class Reader( object ):
    """
    Iterate over all maf blocks in a file in order
    """
    def __init__( self, file, **kwargs ):
        self.file = file
        self.maf_kwargs = kwargs
        # Read and verify maf header, store any attributes
        fields = self.file.readline().split()
        if fields[0] != '##maf': raise "File does not have MAF header"
        self.attributes = parse_attributes( fields[1:] )

    def next( self ):
        return read_next_maf( self.file, **self.maf_kwargs )

    def __iter__( self ):
        return ReaderIter( self )

    def close( self ):
        self.file.close()

class ReaderIter( object ):
    """
    Adapts a `Reader` to the iterator protocol.
    """
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
            # "Empty component" generates an 'e' row
            if c.empty:
                rows.append( ( "e", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.synteny_empty ) )
                continue
            # Regular component
            rows.append( ( "s", c.src, str( c.start ), str( c.size ), c.strand, str( c.src_size ), c.text ) )
            # If component has quality, write a q row
            if c.quality is not None:
                rows.append( ( "q", c.src, "", "", "", "", c.quality ) )
            # If component has synteny follow up with an 'i' row
            if c.synteny_left and c.synteny_right:
                rows.append( ( "i", c.src, "", "", "", "", " ".join( map( str, c.synteny_left + c.synteny_right ) ) ) )
        self.file.write( format_tabular( rows, "llrrrrl" ) )
        self.file.write( "\n" )

    def close( self ):
        self.file.close()

# ---- Helper methods -------------------------------------------------------

def from_string( string, **kwargs ):
    return read_next_maf( StringIO( string ), **kwargs )

def read_next_maf( file, species_to_lengths=None, parse_e_rows=False ):
    """
    Read the next MAF block from `file` and return as an `Alignment` 
    instance. If `parse_i_rows` is true, empty components will be created 
    when e rows are encountered.
    """
    alignment = Alignment(species_to_lengths=species_to_lengths)
    # Attributes line
    line = readline( file, skip_blank=True )
    if not line: return None
    fields = line.split() 
    if fields[0] != 'a': raise "Expected 'a ...' line"
    alignment.attributes = parse_attributes( fields[1:] )
    if 'score' in alignment.attributes:
        alignment.score = alignment.attributes['score']
        del alignment.attributes['score']
    else:
        alignment.score = 0
    # Sequence lines
    last_component = None
    while 1:
        line = readline( file )
        # EOF or Blank line terminates alignment components
        if not line or line.isspace(): break
        if line.isspace(): break 
        # Parse row
        fields = line.split()
        if fields[0] == 's':
            # An 's' row contains sequence for a component
            component = Component()
            component.src = fields[1]
            component.start = int( fields[2] )
            component.size = int( fields[3] )
            component.strand = fields[4]
            component.src_size = int( fields[5] )
            if len(fields) > 6: component.text = fields[6].strip()
            # Add to set
            alignment.add_component( component )
            last_component = component
        elif fields[0] == 'e':
            # An 'e' row, when no bases align for a given species this tells
            # us something about the synteny 
            if parse_e_rows:
                component = Component()
                component.empty = True
                component.src = fields[1]
                component.start = int( fields[2] )
                component.size = int( fields[3] )
                component.strand = fields[4]
                component.src_size = int( fields[5] )
                component.text = None
                synteny = fields[6].strip()
                assert len( synteny ) == 1, \
                    "Synteny status in 'e' rows should be denoted with a single character code"
                component.synteny_empty = synteny
                alignment.add_component( component )
                last_component = component
        elif fields[0] == 'i':
            # An 'i' row, indicates left and right synteny status for the 
            # previous component, we hope ;)
            assert fields[1] == last_component.src, "'i' row does not follow matching 's' row"
            last_component.synteny_left = ( fields[2], int( fields[3] ) )
            last_component.synteny_right = ( fields[4], int( fields[5] ) )
        elif fields[0] == 'q':
            assert fields[1] == last_component.src, "'q' row does not follow matching 's' row"
            # TODO: Should convert this to an integer array?
            last_component.quality = fields[2]
            
    return alignment

def readline( file, skip_blank=False ):
    """Read a line from provided file, skipping any blank or comment lines"""
    while 1:
        line = file.readline()
        #print "every line: %r" % line
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
        
