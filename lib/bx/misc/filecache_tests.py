import filecache
import os
import random
import sys

"""
T="/Users/james/cache/hg18/align/multiz28way/chr10.maf"

def test():
    s = os.stat( T ).st_size
    real_f = open( T )
    f = filecache.FileCache( real_f, s )
    for i in range( 1000 ):
        f.readline()
        
def test_random_seeking():
    s = os.stat( T ).st_size
    raw = open( T )
    f = filecache.FileCache( open( T ), s )
    for i in range( 10000 ):
        seek_to = random.randrange( s )
        
        f.seek( seek_to )
        raw.seek( seek_to )

        l1 = f.readline()
        l2 = raw.readline()
        
        assert l1 == l2
"""
