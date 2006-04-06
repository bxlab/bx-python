import unittest
import os.path
import bx.seq, fasta_tests, nib_tests, qdna_tests

test_fasta = fasta_tests.test_seq
test_nib   = nib_tests.test_seq
test_qdna  = qdna_tests.test_seq

class SEQTestCase (unittest.TestCase):

    def test_get_fasta (self):
        fastafile = bx.seq.seq_file (file (os.path.join('lib','bx','seq','test.fa'),"rb"))
        do_test_get (fastafile, test_fasta, 3, 40)

    def test_get_nib (self):
        nibfile = bx.seq.seq_file (file (os.path.join('lib','bx','seq','test.nib'),"rb"))
        do_test_get (nibfile, test_nib, 3, 40)

    def test_get_qdna (self):
        qdnafile = bx.seq.seq_file (file (os.path.join('lib','bx','seq','test.qdna'),"rb"))
        do_test_get (qdnafile, test_qdna, 3, 40)

def do_test_get (seqfile, test_seq, start, len):
    assert seqfile.get (start, len) == test_seq[start:start+len]

test_classes = [SEQTestCase]
suite = unittest.TestSuite ([unittest.makeSuite (c) for c in test_classes])
