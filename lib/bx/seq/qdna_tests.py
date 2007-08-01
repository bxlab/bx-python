"""
Tests for `bx.seq.qdna`.
"""

import unittest
import sys
import os.path
import qdna

test_qdna = "test_data/seq_tests/test.qdna"

# Same sequence data as stored in test.qdna

valid_seq = "C7wMwHQrMKqEtSREuUv5nsLinpTS8l7jXpbI7IipvCbHnhOdgx" \
          + "5tzRgzYl4j85d:xSlvKPEKEIvZkfiX1YPkBi1Ibhfn9fTZd8gG" \
          + "Wy284hJnwf93W4eHOjeRk7LuVYmH{UTYkYM:b4J4MruMq1ihhv" \
          + "1Yl5W[xXEmi8[JuuLRgooBpy23PllMuUiIiKVIK5vzhjPPYp5Y" \
          + "1eqPxo[e5I24KeCdTV94MZWNybUb:McC:1n4Jczk8JqnR4q1gY" \
          + "HjLS4Bes3s5YvvWdKzS4VrFZy2erhd7YoWRoS[UK8JtSp1{Z1o" \
          + "5:TpvN8mrmWrghiNw{S6nT8DSfF{1ff6kNGpI:FsZE2RgipTVO" \
          + "mJN6vPm8MUgNYd7MDBEu37YOPzPjO1dr"

valid_seq_len = len(valid_seq)

class QDNATestCase(unittest.TestCase):

    def test_get(self):
        qdnafile = qdna.QdnaFile(file(test_qdna,"rb"))
        check_get(qdnafile, 0, valid_seq_len)
        check_get(qdnafile, 0, 40)
        check_get(qdnafile, valid_seq_len - 40, 40)

def check_get(qdnafile, start, len):
    assert qdnafile.get(start, len) == valid_seq[start:start+len]

test_classes = [ QDNATestCase ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
