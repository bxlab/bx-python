"""
Tests for `bx.seq.seq`.
"""


import unittest
import os.path
import sys
import bx.seq, fasta_tests, nib_tests, qdna_tests

test_fa     = "test_data/seq_tests/test.fa"
test2_fa    = "test_data/seq_tests/test2.fa"
test_nib    = "test_data/seq_tests/test.nib"
test_qdna   = "test_data/seq_tests/test.qdna"

valid_fasta = fasta_tests.valid_seq
valid_nib   = nib_tests.valid_seq
valid_qdna  = qdna_tests.valid_seq

# Same sequences as stored in test2.fa

valid2_fa = [("apple",      "GGCGCTGCGATAAGGTTGCGACAACACGGACCTTCTTTTGCCTACCTCTGTTCTTGGCACG"),
             ("orange",     "CGTGCCGAGAACAGAAAATACGCCGGGCGGTGCAGTAGTATCTTGGTATCCGATATGCAGG"),
             ("grapefruit", "CCTGCATATCGACTAGTACACCCTCCCGAGGTACCCCACCCATCCCTCTTTTCTCGGCGCG")]

class SEQTestCase (unittest.TestCase):

    def test_get_fasta (self):
        fastafile = bx.seq.seq_file (file (test_fa, "rb"))
        check_get (fastafile, valid_fasta, 3, 40)

    def test_get_nib (self):
        nibfile = bx.seq.seq_file (file (test_nib, "rb"))
        check_get (nibfile, valid_nib, 3, 40)

    def test_get_qdna (self):
        qdnafile = bx.seq.seq_file (file (test_qdna, "rb"))
        check_get (qdnafile, valid_qdna, 3, 40)

    def test_get_reader (self):
        reader = bx.seq.seq_reader (file (test2_fa, "rb"))
        for (ix,seq) in enumerate(reader):
            assert (ix < len(valid2_fa)), "FastaReader returns too many sequences"
            text = "%s" % seq
            fields = text.split()
            assert (len(fields) == 2), "SeqReader.__str__ returns incorrect sequence string \"%s\" (%d)" % text
            assert (fields[0] == valid2_fa[ix][0]), "FastaReader returned the wrong name (%s,%s)" % (fields[0],valid2_fa[ix][0])
            assert (fields[1] == valid2_fa[ix][1]), "FastaReader returned the wrong text (%s,%s)" % (fields[1],valid2_fa[ix][1])

def check_get (seqfile, valid_seq, start, len):
    assert seqfile.get (start, len) == valid_seq[start:start+len]

test_classes = [SEQTestCase]
suite = unittest.TestSuite ([unittest.makeSuite (c) for c in test_classes])
