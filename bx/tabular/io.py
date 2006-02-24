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
    def __str__( self ):
        return "\t".join( self.fields )

class Header( object ):
    """
    Header of a table -- contains column names and a mapping from them
    to column indexes
    """
    def __init__( self, fields ):
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

class Reader( object ):
    """
    Reader for iterating tabular data
    """
    def __init__( self, input, return_header=True, return_comments=True ):
        self.input = input
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
            return TableRow( self, fields )
        except ParseError, e:
            raise ParseError( str( e ) + "on line " + str( self.linenum ) )    
                
suite = doctest.DocTestSuite( sys.modules[ __name__ ] )
        
        
        
        
        
        
        
        
        
        
        
