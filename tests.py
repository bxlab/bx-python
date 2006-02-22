import unittest
import bx.bitset_tests
import bx.align.score_tests
import bx.intervals.io

suite = unittest.TestSuite( [ bx.bitset_tests.suite, 
                              bx.align.score_tests.suite,
                              bx.intervals.io.suite ] )