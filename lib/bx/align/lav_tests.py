"""
Tests for `bx.align.lav`.
"""

import unittest

import bx.align.lav as lav

test_lav = "test_data/lav_tests/apple_orange.lav"


class lavTestCase(unittest.TestCase):
    def testReader(self):
        reader = lav.Reader(open(test_lav))

        a = next(reader)
        assert a.score == 10286, "a.score is wrong: %s" % a.score
        assert len(a.components) == 2
        check_component(
            a.components[0],
            "apple",
            106,
            252,
            "+",
            411,
            "GTCCGGCCGGCTGAGAGCTACAATACACATGCACGCAGTTTGGCCACTCACATTAAGTATATGAGGAAGGGTTAGCATGAGTTGTACTATAAGGCAGCGGATAGCAGGTTGTGGAAAAATATCCTCCCGATTCAAATCCCCAGGTGCCTAAA----------------GTAGGGCCGGTAGTTGAATGCTTGCCTGTCAGACTGGATGACCAAGTTCAGTATCAACACAATATAGTGCCAGGAGCTAATTGTTCCCCAGCAGCGTGAC",
        )
        check_component(
            a.components[1],
            "lav_tests.orange",
            53,
            252,
            "+",
            361,
            "GTCCGGCCGGCTGTGTGCTACAATACACGTTCACGCAGTTTGGCCAATCACTTTAAGTATATACGAAATGGTTACCATGAGTTGTACTGTAAGGCAGCGGAAAGC---TTGTTAA--------CTCCTGGGCGACATT----GGGGCTGCAACATCGTTTATCCTCCTCTACAACCAATAGCTG-TTGCTTCTTGGTTCAAGTATATCCCATGGATTAGTATCAACACGATATAGTGTCAGGAGCTAATTGTTCCCCAGCAGCGTGAC",
        )

        a = next(reader)
        assert a.score == 3586, "a.score is wrong: %s" % a.score
        assert len(a.components) == 2
        check_component(
            a.components[0],
            "apple",
            52,
            72,
            "+",
            411,
            "TGCATATCGACTATTACAGCCACGCGAGTTACATTCCTCTTTTTTTTTGCTGGCGTCCGGCCGGCTGAGAGC",
        )
        check_component(
            a.components[1],
            "lav_tests.orange",
            2,
            72,
            "-",
            361,
            "TGCATATCGACTAGTACAGCCTCTCGAGTTACCCCCCCCATTCCTCTTGCTGACGTCACGCTGCTGGGGAAC",
        )

        a = next(reader)
        assert a is None

        reader.close()


def check_component(c, src, start, size, strand, src_size, text):
    # ..print "\"%s\" == \"%s\"" % (c.src,src)
    assert c.src == src, f"c.src = {c.src} (expected {src})"
    assert c.start == start, f"c.start = {c.start} (expected {start})"
    assert c.size == size, f"c.size = {c.size} (expected {size})"
    assert c.strand == strand, f"c.strand = {c.strand} (expected {strand})"
    assert c.src_size == src_size, f"c.src_size = {c.src_size} (expected {src_size})"
    assert c.text == text, f'c.text = "{c.text}" (expected "{text}")'
