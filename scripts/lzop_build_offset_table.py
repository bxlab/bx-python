#!/usr/bin/env python2.4

"""
Read a compressed file as created by 'lzop' from stdin and write a table to 
stdout containing the blocksize and the start offset (in bytes) of each 
compressed block. 

usage: %prog < FILENAME.lzo > FILENAME.lzot
"""

import struct, sys

MAGIC="\x89\x4c\x5a\x4f\x00\x0d\x0a\x1a\x0a"

F_ADLER32_D     = 0x00000001L
F_ADLER32_C     = 0x00000002L
F_H_EXTRA_FIELD = 0x00000040L
F_H_GMTDIFF     = 0x00000080L
F_CRC32_D       = 0x00000100L
F_CRC32_C       = 0x00000200L
F_MULTIPART     = 0x00000400L
F_H_FILTER      = 0x00000800L
F_H_CRC32       = 0x00001000L

assert struct.calcsize( "!H" ) == 2
assert struct.calcsize( "!I" ) == 4

class UnpackWrapper( object ):
    def __init__( self, file ):
        self.file = file
    def read( self, amt ):
        return self.file.read( amt )
    def get( self, fmt ):
        t = struct.unpack( fmt, self.file.read( struct.calcsize( fmt ) ) )
        return t[0]

def main():
    f = UnpackWrapper( sys.stdin )
    # Read header
    magic = f.read(9)
    assert magic == MAGIC, "Not LZOP file"
    version = f.get( "!H" )
    lib_version = f.get( "!H" )
    if version >= 0x0940:
        extract_version = f.get( "!H" )
    method = f.get( "!B" )
    assert 1 <= method <= 3, "Only LZO compression is currently supported"
    level = f.get( "!B" )
    flags = f.get( "!I" )
    assert not( flags & F_H_FILTER ), "LZOP filters not supported"
    has_compressed_crc = ( flags & F_CRC32_C or flags & F_ADLER32_C )
    has_uncompressed_crc = ( flags & F_CRC32_D or flags & F_ADLER32_D )
    mode = f.get( "!I" )
    time = f.get( "!I" )
    time_offset = f.get( "!I" )
    fname_len = f.get( "!B" )
    fname = f.read( fname_len )
    assert len( fname ) == fname_len, "EOF reading filename"
    header_crc = f.get( "!I" )
    if ( flags & F_H_EXTRA_FIELD ):
        extra_len = f.get( "!I" )
        extra = f.read( extra_len )
        assert len( extra ) == extra_len, "EOF reading extra field"
    # Done with header
    block_size = None
    expect_no_more = False
    # Read blocks
    while 1:
        size = f.get( "!I" )
        if size == 0: break
        assert not( expect_no_more ), \
            "Encountered an undersized block that was not the last block"
        if block_size is None:
            print "s", size
            block_size = size
        else:
            if size < block_size:
                expect_no_more = True
        compressed_size = f.get( "!I" )
        if has_uncompressed_crc:
            crc = f.get( "!I" )
        if has_compressed_crc:
            compressed_crc = f.get( "!I" )
        print "o", f.file.tell(), compressed_size, size
        compressed_data = f.read( compressed_size )
        assert len( compressed_data ) == compressed_size, \
            "EOF reading compressed data"

if __name__ == "__main__":
    main()
