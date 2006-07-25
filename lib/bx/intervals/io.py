import sys
from itertools import *
from bx.tabular.io import *
from bx.bitset import *

class GenomicInterval( TableRow ):
    """
    A genomic interval stored in a set of fields (a row of a table)
    """
    def __init__( self, reader, fields, chrom_col, start_col, end_col, strand_col, default_strand ):
        TableRow.__init__( self, reader, fields )
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
        if strand_col >= nfields or strand < 0:
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
                  default_strand="+", return_header=True, return_comments=True, force_header=None ):
        TableReader.__init__( self, input, return_header, return_comments, force_header )
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.default_strand = default_strand
    def parse_row( self, line ):
        return GenomicInterval( self, line.split( "\t" ), self.chrom_col, 
                                self.start_col, self.end_col,
                                self.strand_col, self.default_strand )

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
