import bx.misc.bgzf

def test_bgzf():
    f = bx.misc.bgzf.BGZFFile( "../test_data/bgzf_tests/test.txt.gz" )
    print f.read( 10 )
    print f.seek( 0 )
    print f.read( 10 )
    
test_bgzf()