import sys
import doctest
from itertools import *

class ParseError( Exception ):
    pass

class GenomicInterval( object ):
    """
    A genomic interval stored in a set of fields (a row of a table)
    """
    def __init__( self, fields, chrom_col, start_col, end_col, strand_col, default_strand ):
        self.value = self.fields = fields
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.nfields = nfields = len( fields )
        # Parse chrom/source column
        if chrom_col >= nfields:
            raise ParseError( "No field for chrom_col (%d)" % chrom_col )
        self.chrom = fields[chrom_col]
        # Parse start column and ensure it is an integer
        if start_col >= nfields:
            raise ParseError( "No field for start_col (%d)" % start_col )
        try:
            self.start = int( fields[start_col] )
        except ValueError, e:
            raise ParseError( "Could not parse start_col: " + str( e ) )
        # Parse end column and ensure it is an integer
        if end_col >= nfields:
            raise ParseError( "No field for end_col (%d)" % end_col )
        try:
            self.end = int( fields[end_col] )
        except ValueError, e:
            raise ParseError( "Could not parse end_col: " + str( e ) )
        # Parse strand and ensure it is valid
        if strand_col >= nfields:
            # This should probable be immutable since the fields are 
            # not updated when it is set
            self.strand = default_strand
        else:
            strand = fields[strand_col]
            if strand not in ( "+", "-"):
                raise ParseError( "Strand must be either '+' or '-'" )
            self.strand = strand
    def __setattr__( self, name, value ):
        if name == "chrom":
            self.fields[self.chrom_col] = str( value )
        elif name == "start":
            self.fields[self.start_col] = str( value )
        elif name == "end":
            self.fields[self.end_col] = str( value )
        elif name == "strand":
            if self.strand_col < self.nfields:
                self.fields[self.strand_col] = str( value )
        object.__setattr__( self, name, value )
    def __str__( self ):
        return "\t".join( self.fields )

class Header( object ):
    """
    Header of a table -- contains column names and a mapping from them
    to column indexes
    """
    def __init__( self, fields ):
        self.fields = fields
        self.fields_to_column = dict( zip( fields, count() ) )
    def __str__( self ):
        return "#" + "\t".join( self.fields )
        
class Comment( object ):
    def __init__( self, line ):
        self.line = line
    def __str__( self ):
        return "#" + self.line

class Reader( object ):
    """
    Reader for iterating a set of intervals in a tab separated file. Can
    also parse header and comment lines if requested.
    
    >>> r = Reader( [ "#chrom\\tname\\tstart\\tend\\textra",
    ...               "chr1\\tfoo\\t1\\t100\\txxx",
    ...               "chr2\\tbar\\t20\\t300\\txxx",
    ...               "#I am a comment",
    ...               "chr2\\tbar\\t20\\t300\\txxx" ], start_col=2, end_col=3 )
    >>> elements = list( r )
    >>> assert type( elements[0] ) is Header
    >>> str( elements[0] )
    '#chrom\\tname\\tstart\\tend\\textra'
    >>> assert type( elements[1] ) is GenomicInterval
    >>> print elements[1].start, elements[1].end
    1 100
    >>> str( elements[1] )
    'chr1\\tfoo\\t1\\t100\\txxx'
    >>> elements[1].start = 30
    >>> print elements[1].start, elements[1].end
    30 100
    >>> str( elements[1] )
    'chr1\\tfoo\\t30\\t100\\txxx'
    >>> assert type( elements[2] ) is GenomicInterval
    >>> assert type( elements[3] ) is Comment
    >>> assert type( elements[4] ) is GenomicInterval
    """
    def __init__( self, input, chrom_col=0, start_col=1, end_col=2, strand_col=5, 
                  default_strand="+", return_header=True, return_comments=True ):
        self.input = input
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.default_strand = default_strand
        self.return_comments = return_comments
        self.return_header = return_header
        self.input_iter = iter( input )
        self.linenum = -1
        self.header = None
    def __iter__( self ):
        return self
    def next( self ):
        line = self.input_iter.next()
        self.linenum += 1
        line = line.rstrip( "\r\n" )
        # Is it a comment line?
        if line.startswith( "#" ):
            line = line[1:]
            # If a comment and the first line we assume it is a header
            if self.linenum == 0:
                fields = line.split( "\t" )
                self.header = Header( fields )
                if self.return_header:
                    return self.header
                else:
                    return self.next()
            else:
                if self.return_comments:
                    return Comment( line )
                else:
                    return self.next()
        # Not a comment, must be an interval
        fields = line.split( "\t" )
        try:
            return GenomicInterval( fields, self.chrom_col, self.start_col, self.end_col, self.strand_col, self.default_strand )
        except ParseError, e:
            raise ParseError( str( e ) + "on line " + str( self.linenum ) )    
                
suite = doctest.DocTestSuite( sys.modules[ __name__ ] )
        
        
        
        
        
        
        
        
        
        
        
        