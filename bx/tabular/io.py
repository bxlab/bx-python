import sys
import doctest
from itertools import *
from UserDict import DictMixin

class ParseError( Exception ):
    pass

class TableRow( object ):
    """
    A row of a table
    """
    def __init__( self, reader, fields ):
        self.reader = reader
        self.fields = fields
    def __getitem__( self, key ):
        if type( key ) == int:
            return self.fields[key]
        elif type( key ) == str:
            if self.reader.header:
                return self.fields[ self.reader.header.field_to_column[key] ]
            else:
                raise TypeError( "column names only supported for files with headers" )
        else:
            raise TypeError( "field indices must be integers or strings" )
    @property
    def fieldnames( self ):
        return self.reader.header.fields
    def __str__( self ):
        return "\t".join( self.fields )

class Header( object ):
    """
    Header of a table -- contains column names and a mapping from them
    to column indexes
    """
    def __init__( self, fields ):
        self.set_fields( fields )
    def set_fields( self, fields ):
        self.fields = fields
        self.field_to_column = dict( zip( fields, count() ) )
    def __getitem__( self, key ):
        if type( key ) == int:
            return self.fields[key]
        elif type( key ) == str:
            if key in self.field_to_column:
                return key
        else:
            raise TypeError( "field indices must be integers or strings" )
    def __str__( self ):
        return "#" + "\t".join( self.fields )
        
class Comment( object ):
    def __init__( self, line ):
        self.line = line
    def __str__( self ):
        return "#" + self.line

class TableReader( object ):
    """
    Reader for iterating tabular data
    """
    def __init__( self, input, return_header=True, return_comments=True, force_header=None ):
        self.input = input
        self.return_comments = return_comments
        self.return_header = return_header
        self.input_iter = iter( input )
        self.linenum = -1
        self.header = force_header
    def __iter__( self ):
        return self
    def next( self ):
        line = self.input_iter.next()
        self.linenum += 1
        line = line.rstrip( "\r\n" )
        # Is it a comment line?
        if line.startswith( "#" ):
            # If a comment and the first line we assume it is a header
            if self.header is None and self.linenum == 0:
                self.header = self.parse_header( line )
                if self.return_header:
                    return self.header
                else:
                    return self.next()
            else:
                if self.return_comments:
                    return self.parse_comment( line )
                else:
                    return self.next()
        # Not a comment, must be an interval
        try:
            return self.parse_row( line )
        except ParseError, e:
            raise ParseError( str( e ) + "on line " + str( self.linenum ) ) 
    def parse_header( self, line ):
        fields = line[1:].split( "\t" )
        return Header( fields )
    def parse_comment( self, line ):
        return Comment( line[1:] )
    def parse_row( self, line ):
        return TableRow( line.split( "\t" ) )
                
suite = doctest.DocTestSuite( sys.modules[ __name__ ] )
        
        
        
        
        
        
        
        
        
        
        
