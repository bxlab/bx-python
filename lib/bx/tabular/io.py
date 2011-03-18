"""
Reading and writing delimited data files (with headers and comments).
"""

import sys
from itertools import *
from UserDict import DictMixin

FIRST_LINE_IS_HEADER = object()

class ParseError( Exception ):
    def __init__( self, *args, **kwargs ):
        Exception.__init__( self, *args )
        self.linenum = kwargs.get("linenum",None)
    def __str__( self ):
        if self.linenum:
            return Exception.__str__( self ) + " on line " + str(self.linenum)
        else:
            return Exception.__str__( self )

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
        if self.line.startswith("#"):
            return self.line
        return "#" + self.line

class TableReader( object ):
    """
    Reader for iterating tabular data
    """
    def __init__( self, input, return_header=True, return_comments=True, force_header=None, comment_lines_startswith = ["#"] ):
        self.input = input
        self.return_comments = return_comments
        self.return_header = return_header
        self.input_iter = iter( input )
        self.linenum = 0
        self.header = force_header
        self.comment_lines_startswith = comment_lines_startswith
    def __iter__( self ):
        return self
    def next( self ):
        line = self.input_iter.next()
        self.linenum += 1
        line = line.rstrip( "\r\n" )
        # Catch blank lines (throw a warning?)
        # This will end up adding a '#' at the beginning of blank lines
        if line == '':
            if self.return_comments:
                return Comment( line )
            else:
                return self.next()
        # Force header?
        if self.header is FIRST_LINE_IS_HEADER and self.linenum == 1:
            self.header = self.parse_header( line )
            if self.return_header:
                return self.header
            else:
                return self.next()
        # Is it a comment line?
        for comment_line_start in self.comment_lines_startswith:
            if line.startswith( comment_line_start ):
                # If a comment and the first line we assume it is a header
                if self.header is None and self.linenum == 1:
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
            e.linenum = self.linenum
            raise e
    def parse_header( self, line ):
        if line.startswith("#"):
            fields = line[1:].split( "\t" )
        else:
            fields = line.split( "\t" )
        return Header( fields )
    def parse_comment( self, line ):
        return Comment( line )
    def parse_row( self, line ):
        return TableRow( self, line.split( "\t" ) )
        
        
        
        
        
        
        
        
        
        
        
