from align import *

"""
Provides a reader for pairwise alignments in AXT format
"""

class Reader:

    def __init__( self, file ):
        self.file = file

    def next( self ):
        # Read block
        line = readline( self.file )
        if not line: return
        fields = line.split()
        seq1 = readline( self.file )
        seq2 = readline( self.file )
        blank = readline( self.file )
        # Build 2 component alignment
        alignment = Alignment()
        component = Component()
        component.src = fields[1]
        component.start = int( fields[2] ) - 1
        end = int( fields[3] )
        component.size = end - component.start 
        component.strand = "+"
        component.text = seq1.strip()
        alignment.add_component( component )
        component = Component()
        component.src = fields[4]
        component.start = int( fields[5] ) - 1
        end = int( fields[6] )
        component.size = end - component.start
        component.strand = fields[7]
        component.text = seq2.strip()
        
        return alignment

    def __iter__( self ):
        while 1:
            v = self.next()
            if not v: break
            yield v

    def close( self ):
        self.file.close()

# Helper methods

def readline( file, skip_blank=False ):
    """Read a line from provided file, skipping any blank or comment lines"""
    while 1:
        line = file.readline()
        if not line: return None 
        if line[0] != '#' and not ( skip_blank and line.isspace() ):
            return line
