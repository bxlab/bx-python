"""
Classes to support "quantum-DNA" files.

:Author: Bob Harris (rsharris@bx.psu.edu)

A quantum DNA sequence is a sequence of bytes, each representing a probability
distribution (vector) over A, C, G, T.  The QdnaFile class encapsulates the
sequence of bytes, while the mapping from byte value to probability vector is
encapsulated by the QdnaCodebook class.

qdna file format
~~~~~~~~~~~~~~~~

Fields can be in big- or little-endian format;  they must match the endianess
of the magic number.

============ ===========   ======================================================
offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
offset 0x04: 00 00 02 00   version 2.0 (fourth byte is sub version)
offset 0x08: 00 00 00 14   header length (in bytes, including this field)
offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
offset 0x18: xx xx xx xx   (for version >= 2.0) P, offset to named
                           .. properties, 0 indicates no properties
offset    N: ...           name (zero-terminated string)
offset    S: ...           data sequence
offset    P: ...           named properties (see below)
============ ===========   ======================================================

The named properties section consists of a list of pairs of zero-terminated
strings.  The list itself is terminated by an empty string (i.e. another
zero).  In each pair, the first is the name of the property and the second
is its value.  Some names are recognized and handled in some specific manner
(see list below this paragraph).  Any unrecognized name is simply added as
an instance variable with that name, as long as it is not already an instance
variable (in which case it is an error).

Recognized properties (at present only one):
  - codebook: A string in qdna code file format (see QdnaCodebook class for details).
"""

from bx.seq.seq import SeqFile,SeqReader
import sys, struct, string
from StringIO import StringIO

qdnaMagic     = 0xC4B47197L    # big endian magic number for qdna files
qdnaMagicSwap = 0x9771B4C4L

class QdnaFile(SeqFile):

    def __init__(self, file, revcomp=False, name="", gap=None, codebook=None):
        SeqFile.__init__(self,file,revcomp,name,gap)
        if (gap == None): self.gap = chr(0)
        assert (revcomp == False), "reverse complement is not supported for qdna files"
        self.codebook = codebook

        self.byte_order = ">"
        magic = struct.unpack(">L", file.read(4))[0]
        if (magic != qdnaMagic):
            if (magic == qdnaMagicSwap):
                self.byte_order = "<"
            else:
                raise "not a quantum-dna file (magic=%08X)" % magic

        self.magic = magic

        # process header

        self.version = struct.unpack("%sL" % self.byte_order,
                                     self.file.read(4))[0]
        if (self.version not in [0x100,0x200]):
            raise "unsupported quantum-dna (version=%08X)" % self.version

        self.headerLength = struct.unpack("%sL" % self.byte_order,
                                          self.file.read(4))[0]
        if (self.headerLength < 0x10):
            raise "unsupported quantum-dna (header len=%08X)" % self.headerLength
        if (self.version == 0x100) and (self.headerLength != 0x10):
            raise "unsupported quantum-dna (version 1.0 header len=%08X)" % self.headerLength

        self.seqOffset  = struct.unpack("%sL" % self.byte_order,
                                        self.file.read(4))[0]
        self.nameOffset = struct.unpack("%sL" % self.byte_order,
                                        self.file.read(4))[0]
        self.length     = struct.unpack("%sL" % self.byte_order,
                                        self.file.read(4))[0]

        self.propOffset = 0
        if (self.headerLength >= 0x14):
            self.propOffset = struct.unpack("%sL" % self.byte_order,
                                            self.file.read(4))[0]

        self.name = ""
        if (self.nameOffset != 0):
            self.file.seek(self.nameOffset)
            self.name = self.read_string()

        if (self.propOffset != 0):
            self.file.seek(self.propOffset)
            while (True):
                name  = self.read_string()
                if (len(name) == 0): break
                value = self.read_string()
                self.set_property(name,value)


    def set_property(self,name,value):
        if (name == "codebook"):
            self.codebook = QdnaCodebook(StringIO(value))
        else:
            raise "named properties as instance variables are not implemented yet"
            # $$$ do this by adding a properties dict and __getitem__/__setitem__
            # $$$ also need to write properties in QdnaWriter.write()


    def read_string(self):
        s = ""
        while (True):
            ch = self.file.read(1)
            if (ch == chr(0)): break
            s += ch
        return s


    def raw_fetch(self, start, length):
        self.file.seek(self.seqOffset + start)
        return self.file.read(length)


    def get_quantum(self, start, length):
        assert (self.codebook != None), \
                "qdna sequence %s has no code book" % self.name
        return [self.codebook[codeNum] for codeNum in self.raw_fetch(start,length)]


class QdnaReader(SeqReader):

    def __init__(self, file, revcomp=False, name="", gap=None, codebook=None):
        SeqReader.__init__(self,file,revcomp,name,gap)
        self.codebook = codebook

    def next(self):
        if (self.seqs_read != 0): return  # qdna files have just one sequence
        seq = QdnaFile(self.file,self.revcomp,self.name,self.gap,self.codebook)
        self.seqs_read += 1
        return seq

"""
A QdnaCodebook maps code numbers to the corresponding probability vector.  The
latter is a hash from symbols (usually "A", "C", "G", or "T") to the
corresponsing probability.  Note that code numbers are of type string.

qdna code file format:

   The file is ascii text and looks something like what's shown below.  Lines
   beginning with # are comments, and columns are assumed to represent A, C, G
   and T (in that order).  Anything other than five columns is an error.  Note
   that code number zero is usually reserved for gaps in quantum sequences, and
   thus usually won't appear in a code file.  Note that code numbers are
   two-digit hexadecimal (to match the textual displays of quantum sequences).

      01  0.111002  0.072588  0.127196  0.689214
      02  0.081057  0.023799  0.098657  0.796487
      03  0.000260  0.003823  0.000336  0.995581
       ... more lines, usually a total of 255 ...
      FF  0.465900  0.008602  0.482301  0.043197
"""

class QdnaCodebook(object):

    def __init__(self,file):
        (self.alphabet,self.codeToProbs) = self.read_codebook(file)


    def __str__(self):
        codeSet = [codeNum for codeNum in self.codeToProbs]
        codeSet.sort()
        return "\n".join([self.vector_text(codeNum) for codeNum in codeSet])

    def vector_text(self,codeNum):
        if (codeNum in self.codeToProbs): vec = self.codeToProbs[codeNum]
        else:                             vec = {}
        for sym in self.alphabet:
            if (sym not in vec):
                vec[sym] = 0.0
        return ("%02X\t" % ord(codeNum)) \
            + "\t".join(["%.6f" % vec[sym] for sym in self.alphabet])


    def __getitem__ (self,codeNum):
        return self.codeToProbs[codeNum]


    def __setitem__ (self,codeNum,value):
        self.codeToProbs[codeNum] = value # value should be hash from symbol
                                          # .. to probability


    def read_codebook(self,codeF):
        alphabet = "ACGT"
        codeToProbs = {}

        for (lineNum,line) in enumerate (codeF):
            lineNum += 1
            line = line.rstrip()
            stripped = line.strip()
            if (stripped == "") or (stripped.startswith("#")):
                continue

            fields = line.split(None)
            if (len(fields) != 5):
                raise "wrong vector size (line %d)" % lineNum

            try:
                codeNum = int(fields[0],16)
            except:
                raise "bad character code %s (line %d)" \
                    % (fields[0],lineNum)

            if (not 0 <= codeNum <= 255):
                raise "character code %s is outside the valid range (line %d)" \
                     % (fields[0],lineNum)

            if (chr(codeNum) in codeToProbs):
                raise "character code %s appears more than once (line %d)" \
                     % (fields[0],lineNum)

            try:
                vec = {}
                for ix in range(1,5):
                    p = float(fields[ix])
                    if (p < 0) or (p > 1): raise ValueError
                    vec[alphabet[ix-1]] = p
            except:
                raise "%s is a bad probability value (line %d)" \
                     % (fields[ix],lineNum)

            codeToProbs[chr(codeNum)] = vec

        return (alphabet,codeToProbs)


class QdnaWriter(object):

    def __init__(self,file):
        self.file = file

    def write(self,seq):
        text = seq.text
        if (text == None): text = ""

        version   = 0x200
        headerLen = 0x014
        offset    = headerLen + 8

        nameOffset = 0
        if (seq.name != None) and (seq.name != ""):
            nameOffset =  0x01C
            offset     += len(seq.name) + 1
            name       =  seq.name + chr(0)

        dataOffset =  offset
        offset     += len(text)

        assert (seq.codebook == None), \
               "QdnaWriter.write() does not support codebooks yet"
        propOffset = 0

        self.file.write(struct.pack("%sL" % seq.byte_order,qdnaMagic))
        self.file.write(struct.pack("%sL" % seq.byte_order,version))
        self.file.write(struct.pack("%sL" % seq.byte_order,headerLen))
        self.file.write(struct.pack("%sL" % seq.byte_order,dataOffset))
        self.file.write(struct.pack("%sL" % seq.byte_order,nameOffset))
        self.file.write(struct.pack("%sL" % seq.byte_order,len(text)))
        self.file.write(struct.pack("%sL" % seq.byte_order,propOffset))
        if (nameOffset != 0): self.file.write(name)
        self.file.write(text)


    def close(self):
        self.file.close()

