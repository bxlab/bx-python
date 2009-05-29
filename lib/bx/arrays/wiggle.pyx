"""
Support for scores in the `wiggle`_ file format used by the UCSC Genome 
Browser.

The positions in the wiggle format are 1-relative, however,
the positions returned match the BED/interval format which is zero-based, half-open.

.. _wiggle: http://genome.ucsc.edu/goldenPath/help/wiggle.html
"""

def parse_header( line ):
    return dict( [ field.split( '=' ) for field in line.split()[1:] ] )

cdef enum linemode:
    MODE_BED
    MODE_FIXED
    MODE_VARIABLE

cdef class IntervalReader:
    """
    Iterator yielding chrom, start, end, strand, value.
    Values are zero-based, half-open.
    Regions which lack a score are ignored.
    """
    cdef object file
    cdef object current_chrom
    cdef long current_pos
    cdef long current_step
    cdef long current_span
    cdef linemode mode
    def __init__( self, file ):
        self.file = file
        self.current_chrom = None
        self.current_pos = -1
        self.current_step = -1
        self.current_span = -1
        self.mode = MODE_BED

    def __iter__( self ):
        return self

    def __next__( self ):
        while 1:
            line = self.file.readline()
            if not line:
                raise StopIteration()
            if line.isspace():
                continue    
            if line[0] == "#":
                continue
            if line[0].isalpha():
                if line.startswith( "track" ) or line.startswith( "browser" ):
                    continue
                elif line.startswith( "variableStep" ):
                    header = parse_header( line )
                    self.current_chrom = header['chrom']
                    self.current_pos = -1
                    self.current_step = -1
                    if 'span' in header:
                        current_span = int( header['span'] )
                    else:
                        current_span = 1
                    self.mode = MODE_VARIABLE
                    continue
                elif line.startswith( "fixedStep" ):
                    header = parse_header( line )
                    self.current_chrom = header['chrom']
                    self.current_pos = int( header['start'] ) - 1
                    self.current_step = int( header['step'] )
                    if 'span' in header:
                        self.current_span = int( header['span'] )
                    else:
                        self.current_span = 1
                    self.mode = MODE_FIXED
                    continue
            elif self.mode == MODE_BED:
                fields = line.split()
                if len( fields ) > 3:
                    if len( fields ) > 5:
                        return fields[0], int( fields[1] ), int( fields[2] ), fields[5], float( fields[3] )
                    else:
                        return fields[0], int( fields[1] ), int( fields[2] ), "+", float( fields[3] )
            elif self.mode == MODE_VARIABLE: 
                fields = line.split()
                pos = int( fields[0] ) - 1
                return self.current_chrom, pos, pos + self.current_span, "+", float( fields[1] )
            elif self.mode == MODE_FIXED:
                return self.current_chrom, self.current_pos, self.current_pos + self.current_span, "+", float( line.split()[0] )
                self.current_pos += self.current_step
            else:
                raise "Unexpected input line: %s" % line.strip()


#class Reader( object ):
#    """
#    Iterator yielding chrom, position, value.
#    Values are zero-based.
#    Regions which lack a score are ignored.
#    """
#    def __init__( self, f ):
#        self.file = f
#        
#    def __iter__( self ):
#        for chrom, start, end, strand, val in IntervalReader( self.file ):
#            for pos in xrange( start, end ):
#                yield chrom, pos, val
