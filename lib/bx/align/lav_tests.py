"""
Tests for `bx.align.lav`.
"""

import unittest
import sys
import bx.align as align
import bx.align.lav as lav

test_lav = "test_data/lav_tests/apple_orange.lav"

class lavTestCase(unittest.TestCase):

    def testReader(self):

        reader = lav.Reader(file(test_lav))

        a = reader.next()
        assert a.score == 10286, "a.score is wrong: %s" % a.score
        assert len(a.components) == 2
        check_component(a.components[0], "apple",            106, 252, "+", 411, "GTCCGGCCGGCTGAGAGCTACAATACACATGCACGCAGTTTGGCCACTCACATTAAGTATATGAGGAAGGGTTAGCATGAGTTGTACTATAAGGCAGCGGATAGCAGGTTGTGGAAAAATATCCTCCCGATTCAAATCCCCAGGTGCCTAAA----------------GTAGGGCCGGTAGTTGAATGCTTGCCTGTCAGACTGGATGACCAAGTTCAGTATCAACACAATATAGTGCCAGGAGCTAATTGTTCCCCAGCAGCGTGAC")
        check_component(a.components[1], "lav_tests.orange",  53, 252, "+", 361, "GTCCGGCCGGCTGTGTGCTACAATACACGTTCACGCAGTTTGGCCAATCACTTTAAGTATATACGAAATGGTTACCATGAGTTGTACTGTAAGGCAGCGGAAAGC---TTGTTAA--------CTCCTGGGCGACATT----GGGGCTGCAACATCGTTTATCCTCCTCTACAACCAATAGCTG-TTGCTTCTTGGTTCAAGTATATCCCATGGATTAGTATCAACACGATATAGTGTCAGGAGCTAATTGTTCCCCAGCAGCGTGAC")

        a = reader.next()
        assert a.score == 3586, "a.score is wrong: %s" % a.score
        assert len(a.components) == 2
        check_component(a.components[0], "apple",             52,  72, "+", 411, "TGCATATCGACTATTACAGCCACGCGAGTTACATTCCTCTTTTTTTTTGCTGGCGTCCGGCCGGCTGAGAGC")
        check_component(a.components[1], "lav_tests.orange",   2,  72, "-", 361, "TGCATATCGACTAGTACAGCCTCTCGAGTTACCCCCCCCATTCCTCTTGCTGACGTCACGCTGCTGGGGAAC")

        a = reader.next()
        assert a is None

        reader.close()

def check_component( c, src, start, size, strand, src_size, text ):
    #..print "\"%s\" == \"%s\"" % (c.src,src)
    assert c.src      == src,      "c.src = %s (expected %s)"          % (c.src,     src)
    assert c.start    == start,    "c.start = %s (expected %s)"        % (c.start,   start)
    assert c.size     == size,     "c.size = %s (expected %s)"         % (c.size,    size)
    assert c.strand   == strand,   "c.strand = %s (expected %s)"       % (c.strand,  strand)
    assert c.src_size == src_size, "c.src_size = %s (expected %s)"     % (c.src_size,src_size)
    assert c.text     == text,     "c.text = \"%s\" (expected \"%s\")" % (c.text,    text)

test_classes = [ lavTestCase ]
suite = unittest.TestSuite([ unittest.makeSuite(c) for c in test_classes ])
