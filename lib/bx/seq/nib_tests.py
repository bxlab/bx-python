"""
Tests for `bx.seq.nib`.
"""

import unittest

from bx.seq import nib

test_nib = "test_data/seq_tests/test.nib"

# Same sequence data as stored in test.nib

valid_seq = (
    "TGGAGGCATTTGTGATTCAATAGATGCAGAAAGAAACCTTCCTAGAGCTG"
    + "GCGTTCTCTAACTAAAAGTGGAAAGTTCTGAGGAATGAGGACTGTTATAA"
    + "ATCCCACCCCACACCGCACCTTCTCCAGGGAAGTTTCATGGCCGTGAAGA"
    + "GGACAGAAAGTGAGAACCAAGATggaactgaataaacaagcttcacactg"
    + "ttagtttccccatatgcttaccttcccacagatgccaaccttggaggcct"
    + "aagaggcctagaatattatcctttgtctgatcatttctctacaaatttat"
    + "tgttctttgttaagatgctacataagcccaaattctaaccacccctttga"
    + "gttacccatcatcaagtttctcccatgtg"
)

valid_seq_len = len(valid_seq)


class NIBTestCase(unittest.TestCase):
    def test_get(self):
        nibfile = nib.NibFile(open(test_nib, "rb"))
        # Try all combinations of even / odd boundaries
        check_get(nibfile, 0, 10)
        check_get(nibfile, 1, 10)
        check_get(nibfile, 0, 11)
        check_get(nibfile, 1, 11)
        # Test near end of file also
        check_get(nibfile, valid_seq_len - 10, 10)
        check_get(nibfile, valid_seq_len - 11, 11)
        # Test really short gets
        check_get(nibfile, 0, 0)
        check_get(nibfile, 1, 0)
        check_get(nibfile, 0, 1)
        check_get(nibfile, 1, 1)
        # Test negative length
        self.assertRaises(AssertionError, nibfile.get, 20, -1)


def check_get(nibfile, start, len):
    assert nibfile.get(start, len) == valid_seq[start : start + len]
