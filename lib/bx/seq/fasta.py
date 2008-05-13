"""
Classes to support FASTA files.

:Author: Bob Harris (rsharris@bx.psu.edu)

A FASTA file contains multiple sequences.  Each sequence is usually DNA.

A typical FASTA file::

   >mule
   TAATACCCCGGATATATGTCCTCACATAGTTCGAGGTCGAGAAAAATGAC
   TTCCCACCAAGTGGACTCAGCTCGAGTAAACGCCAACGATACGTCCATTA
   GGTGTGTGCCgcaactagtcggacccgttgtgacggaaacaggtccccgc
   caagtcacacgggcatgtcatggacTCTCGATCGTTCATCGCCTTCTTGG
   GTACCGCAGCCGCAATTAAGCCGTGTCTTCTTCCCCCTTCAAACGGGAAT
   CGTGTCGACTTCTTAGGAGCAGNNNNNNNNNNCTAACTCCAGAG
   >donkey
   TAATACCCCGGATATATGTCTTAACATAGTTCCAGGTCGAGAAGAATGAC
   TTGCCACCAAGTGGACTCAGATTCAGTCAACGCGAACGATAAGTCCATTA
   GGTGTGTACCgcaactagtgggacccgttgtgacggaaacaggtcaccgc
   caagtcacacgtgcatgtcatgtacTCTCGATCGTTTATCGCCTTCTTGG
   GTACCGCAGCCGAAATTAAGCCGTGTCTTCTTCCCACTTCAAACGGGAAT
   CGTGTCGACTTTACAGGAACAGNNNNNNNNNNATAACGCCAGAG
    ... more sequences

Typical use:

	for seq in bx.seq.fasta.FastaReader(sys.stdin):
		print seq.name
		print seq.get(0,seq.length)

"""

from bx.seq.seq import SeqFile,SeqReader
import sys, string

class FastaFile(SeqFile):

    def __init__(self, file, revcomp=False, name="", gap=None, lookahead=None, contig=None):
        SeqFile.__init__(self,file,revcomp,name,gap)
        self.lookahead = lookahead
        if (contig == None): contig = 1
        assert (contig >= 1), "contig %d is not legal" % contig

        # nota bene: certainly not the most efficient or elegant implementation

        currContig = 1
        while (True):
            if (self.lookahead != None):
                (line,self.lookahead) = (self.lookahead,None)
            else:
                line = self.file.readline()
            if (line == ""): break
            if not line:
                break
            if (line.startswith(">")):
                if (self.text != None):
                    if (currContig == contig):
                        self.lookahead = line # (next sequence header)
                        break
                    currContig += 1
                self.name = self.extract_name(line[1:])
                self.text = []
                continue
            line = line.split() # (remove whitespace)
            if (self.text == None): self.text = line # (allows headerless fasta)
            else:                   self.text.extend(line)
        assert (currContig == contig), \
            "contig %d is not legal (file contains only %d)" % (contig,currContig)
        if (self.text != None):
            self.text   = "".join(self.text)
            self.length = len(self.text)

class FastaReader(SeqReader):
    
    def __init__(self, file, revcomp=False, name="", gap=None):
        SeqReader.__init__(self,file,revcomp,name,gap)
        self.lookahead = None

    def next(self):
        seq = FastaFile(self.file,self.revcomp,self.name,self.gap,self.lookahead)
        if (seq.text == None): return
        self.lookahead = seq.lookahead
        self.seqs_read += 1
        return seq


class FastaWriter(object):

    def __init__(self,file,columns=50):
        self.file    = file
        self.columns = columns

    def write(self,seq):
        print >>self.file,">%s" % seq.name
        text = seq.text
        if (self.columns != None) and (self.columns > 0):
            text = "\n".join([text[ix:ix+self.columns] \
                              for ix in range(0,len(text),self.columns)])
        print >>self.file,text

    def close(self):
        assert (self.file != None)
        self.file.close()
        self.file = None
