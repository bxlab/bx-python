"""
Tools for fusing contiguous alignment blocks together.
"""

from itertools import *
from copy import deepcopy

def fuse_list( mafs ):
    """
    Try to fuse a list of blocks by progressively fusing each adjacent pair.
    """
    last = None
    for m in mafs:
        if last is None:
            last = m
        else:
            fused = fuse( last, m )
            if fused:
                last = fused
            else:
                yield last
                last = m
    if last:
        yield last

def fuse( m1, m2 ):
    """
    Attempt to fuse two blocks. If they can be fused returns a new block, 
    otherwise returns None.
      
    Example:
      
    >>> import bx.align.maf
      
    >>> block1 = bx.align.maf.from_string( '''
    ... a score=0.0
    ... s hg18.chr10 52686 44 + 135374737 GTGCTAACTTACTGCTCCACAGAAAACATCAATTCTGCTCATGC
    ... s panTro1.chrUn_random 208115356 44 - 240967748 GTGCTAACTGACTGCTCCAGAGAAAACATCAATTCTGTTCATGT
    ... ''' )
    
    >>> block2 = bx.align.maf.from_string( '''
    ... a score=0.0
    ... s hg18.chr10 52730 69 + 135374737 GCAGGTACAATTCATCAAGAAAGGAATTACAACTTCAGAAATGTGTTCAAAATATATCCATACTTTGAC
    ... s panTro1.chrUn_random 208115400 69 - 240967748 GCAGCTACTATTCATCAAGAAAGGGATTACAACTTCAGAAATGTGTTCAAAGTGTATCCATACTTTGAT
    ... ''' )
    
    >>> fused = fuse( block1, block2 )
    
    >>> print fused
    a score=0.0
    s hg18.chr10 52686 113 + 135374737 GTGCTAACTTACTGCTCCACAGAAAACATCAATTCTGCTCATGCGCAGGTACAATTCATCAAGAAAGGAATTACAACTTCAGAAATGTGTTCAAAATATATCCATACTTTGAC
    s panTro1.chrUn_random 208115356 113 - 240967748 GTGCTAACTGACTGCTCCAGAGAAAACATCAATTCTGTTCATGTGCAGCTACTATTCATCAAGAAAGGGATTACAACTTCAGAAATGTGTTCAAAGTGTATCCATACTTTGAT
    <BLANKLINE>
    """
    # Check if the blocks are adjacent, return none if not.
    if len( m1.components ) != len( m2.components ): return None
    for c1, c2 in izip( m1.components, m2.components ):
        if c1.src != c2.src: return None
        if c1.strand != c2.strand: return None
        if c1.end != c2.start: return None
    # Try to fuse:
    n = deepcopy( m1 )
    for c1, c2 in izip( n.components, m2.components ):
        c1.text += c2.text
        c1.size += c2.size
    n.text_size = len( n.components[0].text )
    return n
    
class FusingAlignmentWriter( object ):
    """
    Wrapper for an alignment Writer which attempts to fuse adjacent blocks
    """
    def __init__( self, maf_writer ):
        self.maf_writer = maf_writer
        self.last = None
    def write( self, m ):
        if not self.last:
            self.last = m
        else:
            fused = fuse( self.last, m )
            if fused:
                self.last = fused
            else:
                self.maf_writer.write( self.last )
                self.last = m
    def close( self ):
        if self.last: self.maf_writer.write( self.last )
        self.maf_writer.close()
