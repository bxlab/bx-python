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

from bx.seq.seq import SeqFile
import sys, string

class FastaFile(SeqFile):

    def __init__(self, file, revcomp=False, name="", gap=None):
        SeqFile.__init__(self,file,revcomp,name,gap)

        for line in self.file:
            if (line.startswith(">")):
                if (self.text != None):
                    break # $$$ (multiple sequences not supported yet)
                self.name = self.extract_name(line[1:])
                self.text = []
                continue
            line = line.split() # (remove whitespace)
            if (self.text == None): self.text = line # (allows headerless fasta)
            else:                   self.text.extend(line)
        self.text   = "".join(self.text)
        self.length = len(self.text)

