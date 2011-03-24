import bx.misc.bgzf

def test_bgzf():
    f = bx.misc.bgzf.BGZFFile( "test_data/bgzf_tests/test.txt.gz" )
    assert f.read( 10 ) == "begin 644 "
    print f.seek( 0 )
    assert f.read( 10 ) == "begin 644 "