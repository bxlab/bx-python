"""
Classes to support nib files.

:Author: James Taylor (james@bx.psu.edu), Bob Harris (rsharris@bx.psu.edu)

A nib sequence is a sequence of DNA, using the 10 character alphabet A,C,G,T,N
(upper and lower case).  The file is packed as 4 bits per character.

nib file format
---------------

Fields can be in big- or little-endian format;  they must match the endianess
of the magic number.

============ =========== ======================================================
offset 0x00: 6B E9 3D 3A big endian magic number (3A 3D E9 6B => little endian)
offset 0x04: xx xx xx xx length of data sequence (counted in characters)
offset 0x08:  ...        data sequence;  most significant nybble in each
                         byte is first in sequence
============ =========== ======================================================
"""

from __future__ import division

from bx.seq.seq import SeqFile,SeqReader
import sys, struct, string, math

import _nib

NIB_MAGIC_NUMBER = 0x6BE93D3A
NIB_MAGIC_NUMBER_SWAP = 0x3A3DE96B
NIB_MAGIC_SIZE = 4
NIB_LENGTH_SIZE = 4

class NibFile(SeqFile):

    def __init__(self, file, revcomp=False, name="", gap=None):
        SeqFile.__init__(self,file,revcomp,name,gap)

        self.byte_order = ">"
        magic = struct.unpack(">L", file.read(NIB_MAGIC_SIZE))[0]
        if (magic != NIB_MAGIC_NUMBER):
            if magic == NIB_MAGIC_NUMBER_SWAP: self.byte_order = "<"
            else: raise Exception("Not a NIB file")
        self.magic = magic
        self.length = struct.unpack("%sL" % self.byte_order, file.read(NIB_LENGTH_SIZE))[0]

    def raw_fetch(self, start, length):
        # Check parameters
        assert start >= 0, "Start must be greater than 0"
        assert length >= 0, "Length must be greater than 0"
        assert start + length <= self.length, "Interval beyond end of sequence"
        # Read block of bytes containing sequence
        block_start = int(math.floor(start / 2))
        block_end = int(math.floor((start + length - 1) / 2))
        block_len = block_end + 1 - block_start
        self.file.seek(NIB_MAGIC_SIZE + NIB_LENGTH_SIZE + block_start)
        raw = self.file.read(block_len)
        # Unpack compressed block into a character string and return
        return _nib.translate_raw_data( raw, start, length  )

class NibReader(SeqReader):
    
    def __init__(self, file, revcomp=False, name="", gap=None):
        SeqReader.__init__(self,file,revcomp,name,gap)

    def next(self):
        if (self.seqs_read != 0): return  # nib files have just one sequence
        seq = NibFile(self.file,self.revcomp,self.name,self.gap)
        self.seqs_read += 1
        return seq


class NibWriter(object):

    def __init__(self,file):
        self.file = file

    def write(self,seq):
        assert (False), "NibWriter.write() is not implemented yet"

    def close(self):
        self.file.close()

