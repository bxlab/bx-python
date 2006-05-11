import unittest
import os.path
import sys
import bx.seq, fasta_tests, nib_tests, qdna_tests

test_fasta = fasta_tests.test_seq
test_nib   = nib_tests.test_seq
test_qdna  = qdna_tests.test_seq

# Same sequences as stored in test2.fa

test2 = [("apple",      "GGCGCTGCGATAAGGTTGCGACAACACGGACCTTCTTTTGCCTACCTCTGTTCTTGGCACG"),
         ("orange",     "CGTGCCGAGAACAGAAAATACGCCGGGCGGTGCAGTAGTATCTTGGTATCCGATATGCAGG"),
         ("grapefruit", "CCTGCATATCGACTAGTACACCCTCCCGAGGTACCCCACCCATCCCTCTTTTCTCGGCGCG")]

class SEQTestCase (unittest.TestCase):

    def setUp(self):
        self.save = sys.stdout # this causes an AttributeError if any of these
        sys.stdout = None      # .. tests inadvertently print something

    def tearDown(self):
        sys.stdout = self.save

    def test_get_fasta (self):
        fastafile = bx.seq.seq_file (file (os.path.join('lib','bx','seq','test.fa'),"rb"))
        do_test_get (fastafile, test_fasta, 3, 40)

    def test_get_nib (self):
        nibfile = bx.seq.seq_file (file (os.path.join('lib','bx','seq','test.nib'),"rb"))
        do_test_get (nibfile, test_nib, 3, 40)

    def test_get_qdna (self):
        qdnafile = bx.seq.seq_file (file (os.path.join('lib','bx','seq','test.qdna'),"rb"))
        do_test_get (qdnafile, test_qdna, 3, 40)

    def test_get_reader (self):
        reader = bx.seq.seq_reader (file (os.path.join('lib','bx','seq','test2.fa'),"rb"))
        for (ix,seq) in enumerate(reader):
            assert (ix < len(test2)), "FastaReader returns too many sequences"
            text = "%s" % seq
            fields = text.split()
            assert (len(fields) == 2), "SeqReader.__str__ returns incorrect sequence string \"%s\" (%d)" % text
            assert (fields[0] == test2[ix][0]), "FastaReader returned the wrong name (%s,%s)" % (fields[0],test2[ix][0])
            assert (fields[1] == test2[ix][1]), "FastaReader returned the wrong text (%s,%s)" % (fields[1],test2[ix][1])

def do_test_get (seqfile, test_seq, start, len):
    assert seqfile.get (start, len) == test_seq[start:start+len]

test_classes = [SEQTestCase]
suite = unittest.TestSuite ([unittest.makeSuite (c) for c in test_classes])
