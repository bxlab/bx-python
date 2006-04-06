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

from bx.seq import *

import sys, string

class FastaFile(object):

    def __init__(self, file):
        self.name = ""
        self.text = None

        for line in file:
            if (line.startswith(">")):
                if (self.text != None):
                    break # $$$ (multiple sequences not supported yet)
                self.name = line[1:].strip()
                self.text = ""
                continue
            line = "".join(line.split()) # (remove whitespace)
            if (self.text == None): self.text =  line
            else:                   self.text += line
        self.length = len(self.text)

    def get(self, start, length):
        assert (start >= 0)
        assert (start + length - 1 < self.length)
        return self.text[start:start+length]

