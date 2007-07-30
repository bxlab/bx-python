import seeklzop
import os
import random
import sys

"""
T="/Users/james/cache/hg18/align/multiz28way/chr10.maf"
C="/Users/james/cache/hg18/align/multiz28way/chr10.maf.lzo"

def test():
    f = seeklzop.SeekableLzopFile( C, C + "t", block_cache_size=20 )
    for line in f:
        pass
        
def test_random_seeking():
    s = os.stat( T ).st_size
    raw = open( T )
    f = seeklzop.SeekableLzopFile( C, C + "t", block_cache_size=20 )
    for i in range( 1000 ):
        seek_to = random.randrange( s )
        
        f.seek( seek_to )
        raw.seek( seek_to )

        l1 = f.readline()
        l2 = raw.readline()
        
        assert l1 == l2, "%r != %r" % ( l1, l2 )
        assert raw.tell() == f.tell(), "tells not equal"
"""
