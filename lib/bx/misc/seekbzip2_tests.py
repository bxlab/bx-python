import tempfile
import commands
import os
import random

import seekbzip2
import bz2

F="/Users/james/work/seek-bzip2/test_random.dat.bz2"

import sys

# def test_linear_reading():
#     raw_data = bz2.BZ2File( F ).read()
#     f = seekbzip2.SeekableBzip2File( F, F + "t" )
#     chunk = 11
#     pos = 0
#     for i in range( ( len(raw_data) // chunk ) + 1 ):
#         a = raw_data[pos:pos+chunk]
#         b = f.read( chunk )
#         assert a == b
#         pos += chunk
#         
# def test_random_seeking():
#     raw_data = bz2.BZ2File( F ).read()
#     f = seekbzip2.SeekableBzip2File( F, F + "t" )
#     for i in range( 100 ):
#         seek_to = random.randrange( len( raw_data ) - 100 )
#         chunk = random.randrange( 10, 1000000 )
# 
#         f.seek( seek_to )
#         a = f.read( chunk )
#         b = raw_data[ seek_to : seek_to + chunk ]
#         
#         assert a == b, "'%s' != '%s' on %dth attempt" % ( a.encode("hex"), b.encode("hex"), i )
