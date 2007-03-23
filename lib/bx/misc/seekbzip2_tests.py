import tempfile
import commands
import os
import random

import seekbzip2
import bz2

F="/Users/james/work/seek-bzip2/test_random.dat.bz2"
T="/Users/james/cache/hg18/align/multiz3way/chrY.maf.bz2"

import sys

# def test_linear_reading():
#     raw_data = bz2.BZ2File( F ).read()
#     f = seekbzip2.SeekableBzip2File( F, F + "t" )
#     chunk = 1221
#     pos = 0
#     for i in range( ( len(raw_data) // chunk ) + 1 ):
#         a = raw_data[pos:pos+chunk]
#         b = f.read( chunk )
#         assert a == b
#         pos += chunk
#         assert f.tell() == min( pos, len(raw_data) )
#         
# def test_random_seeking():
#     raw_data = bz2.BZ2File( F ).read()
#     f = seekbzip2.SeekableBzip2File( F, F + "t" )
#     for i in range( 10 ):
#         seek_to = random.randrange( len( raw_data ) - 100 )
#         chunk = random.randrange( 10, 20 )
# 
#         f.seek( seek_to )
#         a = f.read( chunk )
#         b = raw_data[ seek_to : seek_to + chunk ]
#         
#         assert a == b, "'%s' != '%s' on %dth attempt" % ( a.encode("hex"), b.encode("hex"), i )
# 
#         assert f.tell() == min( seek_to + chunk, len(raw_data) )
#         
# def test_text_reading():
#     raw_data = bz2.BZ2File( T ).read()
#     f = seekbzip2.SeekableBzip2File( T, T + "t" )
#     raw_lines = raw_data.split( "\n" )
#     pos = 0
#     for i, line in enumerate( f ):
#         assert line.rstrip( "\n" ) == raw_lines[i], "%d: %r != %r" % ( i, line.rstrip( "\n" ), raw_lines[i] )
#         pos += len( line )
#         ftell = f.tell()
#         assert ftell == pos, "%d != %d" % ( ftell, pos )
#         
#   
# def test_text_reading_2():
#     raw_data = bz2.BZ2File( T ).read()
#     f = seekbzip2.SeekableBzip2File( T, T + "t" )
#     raw_lines = raw_data.split( "\n" )
#     pos = 0
#     i = 0
#     while 1:
#         line = f.readline()
#         if line == "": break
#         assert line.rstrip( "\r\n" ) == raw_lines[i], "%r != %r" % ( line.rstrip( "\r\n" ), raw_lines[i] )
#         pos += len( line )
#         ftell = f.tell()
#         assert ftell == pos, "%d != %d" % ( ftell, pos )  
#         i += 1    
        
        