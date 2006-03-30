from bx.align import *

import itertools
from bx import interval_index_file

class Reader( object ):
    """Iterate over all lav blocks in a file in order"""
    
    def __init__( self, file, species_to_lengths=None ):
        self.file = file
        self.species_to_lengths = species_to_lengths
        self.seq1_name = None
        self.seq1_header = None
        self.seq1_start = None
        self.seq1_end = None
        self.seq2_name = None
        self.seq2_header = None
        self.seq2_start = None
        self.seq2_end = None

    def next( self ):
        while 1:
            line = self.file.readline().rstrip()
            if not line: 
                return None
            if line.startswith( "s {" ):
                self.parse_s_stanza()
            elif line.startswith( "h {" ):
                self.parse_h_stanza()
            elif line.startswith( "a {" ):
                return self.parse_a_stanza()
            else:
                pass
                
    def parse_s_stanza( self ):
        fields = self.file.readline().split()
        self.seq1_name = fields[0].strip( '"' )
        self.seq1_start = int( fields[ 1 ] ) - 1
        self.seq1_end = int( fields[ 2 ] )
        fields = self.file.readline().split()
        self.seq2_name = fields[0].strip( '"' )
        self.seq2_start = int( fields[ 1 ] ) - 1
        self.seq2_end = int( fields[ 2 ] )
        
    def parse_h_stanza( self ):
        line = self.file.readline().strip( '"' )
        if line.startswith( ">" ):
            self.seq1_header = line[1:]
        else:
            raise "Non FASTA lav not supported"
        line = self.file.readline().strip( '"' )
        if line.startswith( ">" ):
            self.seq2_header = line[1:]
        else:
            raise "Non FASTA lav not supported"
            
    def parse_a_stanza( self ):
        # 's' line -- score, 1 field
        line = self.file.readline().split()
        assert fields[0] == "s", "s line expected in a stanza"
        score = float( fields[1] )
        # 'b' line -- begin positions in seqs, 2 fields
        line = self.file.readline().split()
        assert fields[0] == "b", "b line expected in a stanza"
        beg1 = int( fields[1] ) - 1
        beg2 = int( fields[2] ) - 1
        # 'e' line -- end positions in seqs, 2 fields
        line = self.file.readline().split()
        assert fields[0] == "b", "b line expected in a stanza"
        len1 = int( fields[1] ) - beg1
        len2 = int( fields[2] ) - beg2
        # 'l' lines
        while 1:
            line = self.file.readline().split()
            if fields[0] != "l":
                break
            start1  = int( fields[1] ) - 1
            start2  = int( fields[2] ) - 1
            length1 = int( fields[3] ) - start1
            length2 = int( fields[4] ) - start2
            pctId   = float( fields[5] )
            

            

    def __iter__( self ):
        return ReaderIter( self )

    def close( self ):
        self.file.close()

class ReaderIter( object ):
    def __init__( self, reader ):
        self.reader = reader
    def __iter__( self ): 
        return self
    def next( self ):
        v = self.reader.next()
        if not v: raise StopIteration
        return v
