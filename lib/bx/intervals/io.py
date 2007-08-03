"""
Support for reading and writing genomic intervals from delimited text files.
"""

import sys
from itertools import *
from bx.tabular.io import *
from bx.bitset import *

class MissingFieldError( ParseError ):
    pass

class FieldFormatError( ParseError ):
    def __init__( self, *args, **kwargs):
        ParseError.__init__( self, *args, **kwargs )
        self.expected = kwargs.get("expected",None)
    def __str__( self ):
        if self.expected:
            return ParseError.__str__( self ) + ", " + self.expected + " expected"
        else:
            return ParseError.__str__( self )

class StrandFormatError( ParseError ):
    pass

class GenomicInterval( TableRow ):
    """
    A genomic interval stored in a set of fields (a row of a table)
    """
    def __init__( self, reader, fields, chrom_col, start_col, end_col, strand_col, default_strand, fix_strand=False ):
        TableRow.__init__( self, reader, fields )
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.nfields = nfields = len( fields )
        # Parse chrom/source column
        if chrom_col >= nfields:
            raise MissingFieldError( "No field for chrom_col (%d)" % chrom_col )
        self.chrom = fields[chrom_col]
        # Parse start column and ensure it is an integer
        if start_col >= nfields:
            raise MissingFieldError( "No field for start_col (%d)" % start_col )
        try:
            self.start = int( fields[start_col] )
        except ValueError, e:
            raise FieldFormatError( "Could not parse start_col: " + str( e ), expected="integer" )
        # Parse end column and ensure it is an integer
        if end_col >= nfields:
            raise MissingFieldError( "No field for end_col (%d)" % end_col )
        try:
            self.end = int( fields[end_col] )
        except ValueError, e:
            raise FieldFormatError( "Could not parse end_col: " + str( e ), expected="integer" )
        # Ensure start <= end
        if self.end < self.start:
            raise ParseError( "Start is greater than End. Interval length is < 1." )
        # Parse strand and ensure it is valid
        if strand_col >= nfields or strand_col < 0:
            # This should probable be immutable since the fields are 
            # not updated when it is set
            self.strand = default_strand
        else:
            strand = fields[strand_col]
            if strand not in ( "+", "-"):
                if fix_strand:
                    strand = "+"
                else: raise StrandFormatError( "Strand must be either '+' or '-'" )
            self.strand = strand
    def __setattr__( self, name, value ):
        if name == "chrom":
            self.fields[self.chrom_col] = str( value )
        elif name == "start":
            self.fields[self.start_col] = str( value )
        elif name == "end":
            self.fields[self.end_col] = str( value )
        elif name == "strand":
            if self.strand_col < self.nfields and self.strand_col >= 0:
                self.fields[self.strand_col] = str( value )
        object.__setattr__( self, name, value )
    def __str__( self ):
        return "\t".join( self.fields )
    def copy( self ):
        return GenomicInterval(self.reader, list( self.fields ), self.chrom_col, self.start_col, self.end_col, self.strand_col, self.strand)

class GenomicIntervalReader( TableReader ):
    """
    Reader for iterating a set of intervals in a tab separated file. Can
    also parse header and comment lines if requested.
    
    >>> r = GenomicIntervalReader( [ "#chrom\\tname\\tstart\\tend\\textra",
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
                  default_strand="+", return_header=True, return_comments=True, force_header=None, fix_strand=False ):
        TableReader.__init__( self, input, return_header, return_comments, force_header )
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.default_strand = default_strand
        self.fix_strand = fix_strand
    def parse_row( self, line ):
        return GenomicInterval( self, line.split( "\t" ), self.chrom_col, 
                                self.start_col, self.end_col,
                                self.strand_col, self.default_strand, fix_strand=self.fix_strand )

    def binned_bitsets( self , upstream_pad=0, downstream_pad=0, lens={} ):
        last_chrom = None
        last_bitset = None
        bitsets = dict()
        for interval in self:
            if type( interval ) == GenomicInterval:
                chrom = interval[self.chrom_col]
                if chrom != last_chrom:
                    if chrom not in bitsets:
                        if chrom in lens:
                            size = lens[chrom]
                        else:
                            size = MAX
                        bitsets[chrom] = BinnedBitSet( size )
                    last_chrom = chrom
                    last_bitset = bitsets[chrom]
                start = max(int( interval[self.start_col]), 0 )
                end = min(int( interval[self.end_col]), size)
                last_bitset.set_range( start, end-start )
        return bitsets

class NiceReaderWrapper( GenomicIntervalReader ):
    def __init__( self, reader, **kwargs ):
        GenomicIntervalReader.__init__( self, reader, **kwargs )
        self.outstream = kwargs.get("outstream", None)
        self.print_delegate = kwargs.get("print_delegate", None)
        self.skipped = 0
        self.skipped_lines = []
        self.input_wrapper = iter( self.input )
        self.input_iter = self.iterwrapper()
    def __iter__( self ):
        return self
    def next( self ):
        while 1:
            try:
                nextitem = GenomicIntervalReader.next( self )
                return nextitem
            except ParseError, e:
                if self.outstream:
                    if self.print_delegate and hasattr(self.print_delegate,"__call__"):
                        self.print_delegate( self.outstream, e, self )
                self.skipped += 1
                # no reason to stuff an entire bad file into memmory
                if self.skipped < 10:
                    self.skipped_lines.append( (self.linenum, self.current_line) )
    def iterwrapper( self ):
        while 1:
            self.current_line = self.input_wrapper.next()
            yield self.current_line
