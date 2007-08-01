"""
Tests for `bx.seq.fasta`.
"""

import unittest
import sys
import os.path
import fasta

test_fa = "test_data/seq_tests/test.fa"

# Same sequence data as stored in test.fa

valid_seq = "TGGAGGCATTTGTGATTCAATAGATGCAGAAAGAAACCTTCCTAGAGCTG" \
          + "GCGTTCTCTAACTAAAAGTGGAAAGTTCTGAGGAATGAGGACTGTTATAA" \
          + "ATCCCACCCCACACCGCACCTTCTCCAGGGAAGTTTCATGGCCGTGAAGA" \
          + "GGACAGAAAGTGAGAACCAAGATggaactgaataaacaagcttcacactg" \
          + "ttagtttccccatatgcttaccttcccacagatgccaaccttggaggcct" \
          + "aagaggcctagaatattatcctttgtctgatcatttctctacaaatttat" \
          + "tgttctttgttaagatgctacataagcccaaattctaaccacccctttga" \
          + "gttacccatcatcaagtttctcccatgtg"

valid_seq_len = len( valid_seq )

class FASTATestCase( unittest.TestCase ):

    def test_get( self ):
        fastafile = fasta.FastaFile( file(test_fa, "rb" ) )
        check_get(fastafile, 0, valid_seq_len)
        check_get(fastafile, 0, 40)
        check_get(fastafile, valid_seq_len - 40, 40)

def check_get( fastafile, start, len ):
    assert fastafile.get( start, len ) == valid_seq[start:start+len]

test_classes = [ FASTATestCase ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
