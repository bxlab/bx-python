from __future__ import print_function

import bx.misc.bgzf
import sys

def test_bgzf():
    print("OH COMEON")
    sys.stdout.flush()
    f = bx.misc.bgzf.BGZFFile( "test_data/bgzf_tests/test.txt.gz" , "r")
    print("OH COMEON")
    sys.stdout.flush()
    assert f.read( 10 ) == "begin 644 "
    print(f.seek( 0 ))
    assert f.read( 10 ) == "begin 644 "
