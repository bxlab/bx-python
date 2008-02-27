"""
Support for scores in the `wiggle`_ file format used by the UCSC Genome 
Browser.

The positions in the wiggle format are 1-relative, however,
the positions returned match the BED/interval format which is zero-based, half-open.

.. _wiggle: http://genome.ucsc.edu/goldenPath/help/wiggle.html
"""

def parse_header( line ):
    return dict( [ field.split( '=' ) for field in line.split()[1:] ] )

def IntervalReader( f ):
    """
    Iterator yielding chrom, start, end, strand, value.
    Values are zero-based, half-open.
    Regions which lack a score are ignored.
    """
    current_chrom = None
    current_pos = None
    current_step = None

    # always for wiggle data
    strand = '+'

    mode = "bed"

    for line in f:
        if line.isspace() or line.startswith( "track" ) or line.startswith( "#" ) or line.startswith( "browser" ):
            continue
        elif line.startswith( "variableStep" ):
            header = parse_header( line )
            current_chrom = header['chrom']
            current_pos = None
            current_step = None
            if 'span' in header: current_span = int( header['span'] )
            else: current_span = 1
            mode = "variableStep"
        elif line.startswith( "fixedStep" ):
            header = parse_header( line )
            current_chrom = header['chrom']
            current_pos = int( header['start'] ) - 1
            current_step = int( header['step'] )
            if 'span' in header: current_span = int( header['span'] )
            else: current_span = 1
            mode = "fixedStep"
        elif mode == "bed":
            fields = line.split()
            if len( fields ) > 3:
                if len( fields ) > 5:
                    yield fields[0], int( fields[1] ), int( fields[2] ), fields[5], float( fields[3] )
                else:
                    yield fields[0], int( fields[1] ), int( fields[2] ), strand, float( fields[3] )
        elif mode == "variableStep": 
            fields = line.split()
            pos = int( fields[0] ) - 1
            yield current_chrom, pos, pos + current_span, strand, float( fields[1] )
        elif mode == "fixedStep":
            yield current_chrom, current_pos, current_pos + current_span, strand, float( line.split()[0] )
            current_pos += current_step
        else:
            raise "Unexpected input line: %s" % line.strip()


class Reader( object ):
    """
    Iterator yielding chrom, position, value.
    Values are zero-based.
    Regions which lack a score are ignored.
    """
    def __init__( self, f ):
        self.file = f
        
    def __iter__( self ):
        for chrom, start, end, strand, val in IntervalReader( self.file ):
            for pos in xrange( start, end ):
                yield chrom, pos, val
