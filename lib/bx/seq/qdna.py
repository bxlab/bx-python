"""
Classes to support "quantum-DNA" files
--------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

A quantum DNA sequence is a sequence of bytes, each representing a probability
distribution (vector) over A, C, G, T.  The mapping from byte value to
probability vector is external to this module.

qdna file format:

   Fields can be in little- or big-endian format (this is determined from the
   magic number).

   offset 0x00: C4 B4 71 97   new magic number
   offset 0x04: 00 00 01 00   version (second byte will be sub version)
   offset 0x08: 00 00 00 10   header length (in bytes, including this field)
   offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
   offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
   offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
   offset    N:  ...          name (zero-terminated string)
   offset    S:  ...          data sequence
"""

import sys, struct, string

qdnaMagic     = 0xC4B47197L    # big endian magic number for qdna files
qdnaMagicSwap = 0x9771B4C4L

class QdnaFile(object):

    def __init__(self, file):
        self.byte_order = ">" 
        magic = struct.unpack(">L", file.read(4))[0]
        if (magic != qdnaMagic):
            if (magic == qdnaMagicSwap):
                self.byte_order = "<"
            else:
                raise "not a quantum-dna file (magic=%08X)" % magic

        self.file  = file
        self.magic = magic

        # process header

        self.version = struct.unpack("%sL" % self.byte_order, 
                                     self.file.read(4))[0]
        if (self.version != 0x100):
            raise "unsupport quantum-dna (version=%08X)" % self.version

        self.headerLength = struct.unpack("%sL" % self.byte_order, 
                                          self.file.read(4))[0]
        if (self.headerLength < 0x10):
            raise "unsupport quantum-dna (header len=%08X)" % self.headerLength

        self.seqOffset  = struct.unpack("%sL" % self.byte_order, 
                                        self.file.read(4))[0]
        self.nameOffset = struct.unpack("%sL" % self.byte_order, 
                                        self.file.read(4))[0]
        self.length     = struct.unpack("%sL" % self.byte_order, 
                                        self.file.read(4))[0]

        self.name = ""
        if (self.nameOffset != 0):
            self.file.seek(self.nameOffset)
            while (True):
                ch = self.file.read(1)
                if (ch == chr(0)): break
                self.name += ch

    def get(self, start, length):
        assert (start >= 0)
        assert (start + length - 1 < self.length)
        self.file.seek(self.seqOffset + start)
        return self.file.read(length)

