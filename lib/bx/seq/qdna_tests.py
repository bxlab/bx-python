import unittest
import os.path
import qdna

# Same sequence data as stored in test.qdna

test_seq = "C7wMwHQrMKqEtSREuUv5nsLinpTS8l7jXpbI7IipvCbHnhOdgx" \
         + "5tzRgzYl4j85d:xSlvKPEKEIvZkfiX1YPkBi1Ibhfn9fTZd8gG" \
         + "Wy284hJnwf93W4eHOjeRk7LuVYmH{UTYkYM:b4J4MruMq1ihhv" \
         + "1Yl5W[xXEmi8[JuuLRgooBpy23PllMuUiIiKVIK5vzhjPPYp5Y" \
         + "1eqPxo[e5I24KeCdTV94MZWNybUb:McC:1n4Jczk8JqnR4q1gY" \
         + "HjLS4Bes3s5YvvWdKzS4VrFZy2erhd7YoWRoS[UK8JtSp1{Z1o" \
         + "5:TpvN8mrmWrghiNw{S6nT8DSfF{1ff6kNGpI:FsZE2RgipTVO" \
         + "mJN6vPm8MUgNYd7MDBEu37YOPzPjO1dr"

test_seq_len = len(test_seq)

class QDNATestCase(unittest.TestCase):

    def test_get(self):
        qdnafile = qdna.QdnaFile(file(os.path.join('lib','bx','seq','test.qdna')))
        do_test_get(qdnafile, 0, test_seq_len)
        do_test_get(qdnafile, 0, 40)
        do_test_get(qdnafile, test_seq_len - 40, 40)

def do_test_get(qdnafile, start, len):
    assert qdnafile.get(start, len) == test_seq[start:start+len]

test_classes = [ QDNATestCase ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
