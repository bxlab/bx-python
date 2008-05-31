"""
Classes for reading and writing motif data.
"""

from bx.motif.pwm import FrequencyMatrix

class TransfacMotif( object ):
    
    def __init__( self ):
        self.accession = None
        self.id = None
        self.dates = None
        self.name = None
        self.description = None
        self.binding_factors = None
        self.basis = None
        self.comment = None
        self.matrix = None
        
class TransfacReader( object ):
    """
    Reads motifs in TRANSFAC format.
    """
    
    def __init__( self, input ):
        self.input = iter( input )
        self.input_exhausted = False
    
    def as_dict( self, key="id" ):
        """
        Return a dictionary containing all remaining motifs, using `key`
        as the dictionary key.
        """
        rval = {}
        for motif in self:
            rval[ getattr( motif, key ) ] = motif
        return rval
    
    def __iter__( self ):
        return self
    
    def next( self ):
        rval = self.next_motif()
        while rval is None:
            rval = self.next_motif()
        return rval
    
    def next_motif( self ):
        if self.input_exhausted:
            raise StopIteration
        # Accumulate lines until either the end of record indicator "//" is
        # encounted or the input is exhausted.   
        lines = []
        while 1:
            try:
                line = self.input.next()
            except StopIteration, e:
                self.input_exhausted = True
                break
            if line.startswith( "//" ):
                break
            if not line.isspace():
                lines.append( line )
        if lines:
            return self.parse_record( lines )
    
    parse_actions = {
        "AC": ( "store_single", "accession" ),
        "ID": ( "store_single", "id" ),
        "DT": ( "store_single_list", "dates" ),
        "NA": ( "store_single", "name" ),
        "DE": ( "store_block", "description" ),
        "BF": ( "store_single_list", "binding_factors" ),
        "BA": ( "store_block", "basis" ),
        "CC": ( "store_block", "comment" ),
        "P0": ( "store_matrix", "matrix" )
    }
    
    def parse_record( self, lines ):
        """
        Parse a TRANSFAC record out of `lines` and return a motif.
        """
        # Break lines up
        temp_lines = []
        for line in lines:
            fields = line.rstrip( "\r\n" ).split( None, 1 )
            if len( fields ) == 1:
                fields.append( "" )
            temp_lines.append( fields )
        lines = temp_lines
        # Fill in motif from lines
        motif = TransfacMotif()
        current_line = 0
        while 1:
            # Done parsing if no more lines to consume
            if current_line >= len( lines ):
                break
            # Remove prefix and first separator from line
            prefix, rest = lines[ current_line ]
            # No action for this prefix, just ignore the line
            if prefix not in self.parse_actions:
                current_line += 1
                continue
            # Get action for line
            action = self.parse_actions[ prefix ]
            # Store a single line value
            if action[0] == "store_single":
                key = action[1]
                setattr( motif, key, rest )
                current_line += 1
            # Add a single line value to a list
            if action[0] == "store_single_list":
                key = action[1]
                if not getattr( motif, key ):
                    setattr( motif, key, [] )
                getattr( motif, key ).append( rest )
                current_line += 1
            # Store a block of text
            if action[0] == "store_block":
                key = action[1]
                value = []
                while current_line < len( lines ) and lines[ current_line ][0] == prefix:
                    value.append( lines[current_line][1] )
                    current_line += 1
                setattr( motif, key, str.join( "\n", value ) )
            # Store a matrix
            if action[0] == "store_matrix":
                # First line is alphabet
                alphabet = rest.split()
                alphabet_size = len( alphabet )
                rows = []
                pattern = ""
                current_line += 1
                # Next lines are the rows of the matrix (we allow 0 rows)
                while current_line < len( lines ):
                    prefix, rest = lines[ current_line ]
                    # Prefix should be a two digit 0 padded row number
                    if not prefix.isdigit():
                        break
                    # The first `alphabet_size` fields are the row values
                    values = rest.split()
                    rows.append( map( float, values[:alphabet_size] ) )
                    # TRANSFAC includes an extra column with the IUPAC code
                    if len( values ) > alphabet_size:
                        pattern += values[alphabet_size]
                    current_line += 1
                # Only store the pattern if it is the correct length (meaning
                # that every row had an extra field)
                if len( pattern ) != len( rows ):
                    pattern = None
                matrix = FrequencyMatrix.from_rows( alphabet, rows )
                setattr( motif, action[1], matrix )
        # Only return a motif if we saw at least ID or AC or NA
        if motif.id or motif.accession or motif.name:
            return motif