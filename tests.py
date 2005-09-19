import unittest
import bx.bitset_tests
import bx.align.score_tests

suite = unittest.TestSuite( [ bx.bitset_tests.suite, bx.align.score_tests.suite ] )
