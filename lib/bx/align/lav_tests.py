import unittest
import sys
import bx.align as align
import bx.align.lav as lav

from StringIO import StringIO

test_lav = """#:lav
d {
  "blastz.v0(256) lib/bx/align/lav_tests_apple.fa[51,361] lib/bx/align/lav_tests_orange.nib K=2000
     A    C    G    T
    91 -114  -31 -123
  -114  100 -125  -31
   -31 -125  100 -114
  -123  -31 -114   91
   O = 400, E = 30, K = 2000, L = 2000, M = 0"
}
#:lav
s {
  "lib/bx/align/lav_tests_apple.fa" 51 361 0 1
  "lib/bx/align/lav_tests_orange.nib" 1 361 0 1
}
h {
  "> apple"
  "lib/bx/align/lav_tests_orange.nib:1-361"
}
a {
  s 10286
  b 57 54
  e 308 305
  l 57 54 161 158 88
  l 165 159 171 165 71
  l 180 166 194 180 53
  l 199 181 208 190 60
  l 209 207 224 222 56
  l 226 223 308 305 76
}
x {
  n 0
}
#:lav
s {
  "lib/bx/align/lav_tests_apple.fa" 51 361 0 1
  "lib/bx/align/lav_tests_orange.nib-" 1 361 1 1
}
h {
  "> apple"
  "lib/bx/align/lav_tests_orange.nib:1-361 (reverse complement)"
}
a {
  s 3586
  b 3 3
  e 74 74
  l 3 3 74 74 72
}
x {
  n 0
}
m {
  n 0
}
#:eof
"""

class lavTestCase(unittest.TestCase):

    def setUp(self):
        sys.stdout = None # this causes an AttributeError if any of these
                          # .. tests inadvertently print something

    def tearDown(self):
        sys.stdout = sys.__stdout__

    def testReader(self):

        reader = lav.Reader(StringIO(test_lav))

        a = reader.next()
        assert a.score == 10286, "a.score is wrong: %s" % a.score
        assert len(a.components) == 2
        check_component(a.components[0], "align.lav_tests_apple",  106, 252, "+", 411, "GTCCGGCCGGCTGAGAGCTACAATACACATGCACGCAGTTTGGCCACTCACATTAAGTATATGAGGAAGGGTTAGCATGAGTTGTACTATAAGGCAGCGGATAGCAGGTTGTGGAAAAATATCCTCCCGATTCAAATCCCCAGGTGCCTAAA----------------GTAGGGCCGGTAGTTGAATGCTTGCCTGTCAGACTGGATGACCAAGTTCAGTATCAACACAATATAGTGCCAGGAGCTAATTGTTCCCCAGCAGCGTGAC")
        check_component(a.components[1], "align.lav_tests_orange",  53, 252, "+", 361, "GTCCGGCCGGCTGTGTGCTACAATACACGTTCACGCAGTTTGGCCAATCACTTTAAGTATATACGAAATGGTTACCATGAGTTGTACTGTAAGGCAGCGGAAAGC---TTGTTAA--------CTCCTGGGCGACATT----GGGGCTGCAACATCGTTTATCCTCCTCTACAACCAATAGCTG-TTGCTTCTTGGTTCAAGTATATCCCATGGATTAGTATCAACACGATATAGTGTCAGGAGCTAATTGTTCCCCAGCAGCGTGAC")

        a = reader.next()
        assert a.score == 3586, "a.score is wrong: %s" % a.score
        assert len(a.components) == 2
        check_component(a.components[0], "align.lav_tests_apple",   52,  72, "+", 411, "TGCATATCGACTATTACAGCCACGCGAGTTACATTCCTCTTTTTTTTTGCTGGCGTCCGGCCGGCTGAGAGC")
        check_component(a.components[1], "align.lav_tests_orange",   2,  72, "-", 361, "TGCATATCGACTAGTACAGCCTCTCGAGTTACCCCCCCCATTCCTCTTGCTGACGTCACGCTGCTGGGGAAC")

        a = reader.next()
        assert a is None

        reader.close()

def check_component( c, src, start, size, strand, src_size, text ):
    assert c.src      == src,      "c.src = %s (expected %s)"          % (c.src,     src)
    assert c.start    == start,    "c.start = %s (expected %s)"        % (c.start,   start)
    assert c.size     == size,     "c.size = %s (expected %s)"         % (c.size,    size)
    assert c.strand   == strand,   "c.strand = %s (expected %s)"       % (c.strand,  strand)
    assert c.src_size == src_size, "c.src_size = %s (expected %s)"     % (c.src_size,src_size)
    assert c.text     == text,     "c.text = \"%s\" (expected \"%s\")" % (c.text,    text)

test_classes = [ lavTestCase ]
suite = unittest.TestSuite([ unittest.makeSuite(c) for c in test_classes ])
