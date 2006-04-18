"""
Classes to support FASTA files
------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
:Version: $Revision: $

A FASTA file contains multiple sequences.  Each sequence is usually DNA.

WARNING:  THIS MODULE CURRENTLY ONLY READS THE FIRST SEQUENCE

A typical FASTA file:

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
"""

from bx.seq.seq import SeqFile,SeqReader
import sys, string

class FastaFile(SeqFile):

    def __init__(self, file, revcomp=False, name="", gap=None,lookahead=None):
        SeqFile.__init__(self,file,revcomp,name,gap)
        self.lookahead = None

        while (True):
            if (lookahead != None): (line,lookahead) = (lookahead,None)
            else:                    line = self.file.readline()
            if (line == ""): break
            if (line.startswith(">")):
                if (self.text != None):
                    self.lookahead = line # (next sequence header)
                    break
                self.name = self.extract_name(line[1:])
                self.text = []
                continue
            line = line.split() # (remove whitespace)
            if (self.text == None): self.text = line # (allows headerless fasta)
            else:                   self.text.extend(line)
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



