import interval_index_file
from interval_index_file import Indexes
from tempfile import mktemp
import random

def test_offsets():
    assert interval_index_file.offsets_for_max_size( 512*1024*1024  - 1 ) == [ 512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0 ]

def test_interval_index_file():
    ix = Indexes()
    chrs = []
    for i in range( 5 ):
        intervals = []
        name = "seq%d" % i
        max = random.randint( 0, interval_index_file.MAX )
        # print name, "size", max
        for i in range( 500 ):
            start = random.randint( 0, max )
            end = random.randint( 0, max )
            if end < start:
                end, start = start, end
            ix.add( name, start, end, i, max=interval_index_file.MAX )
            intervals.append( ( start, end, i ) )
        chrs.append( intervals )
    fname = mktemp()
    f = open( fname, "w" )
    ix.write( f )
    f.close()
    del ix
    
    ix = Indexes( fname )
    for i in range( 5 ):
        intervals = chrs[i]
        name = "seq%d" % i
        for i in range( 100 ):
            start = random.randint( 0, max )
            end = random.randint( 0, max )
            if end < start:
                end, start = start, end
            query_intervals = set()
            for ( s, e, i ) in intervals:
                if e > start and s < end:
                    query_intervals.add( ( s, e, i ) )
            result = ix.find( name, start, end )
            for inter in result:
                assert inter in query_intervals

def test_zero():
    ix = Indexes()
    ix.add("t.idx", 0, 0, 1, 123)
