import sys; sys.path.append( "lib" )

import unittest
import bx.bitset_tests
import bx.align.score_tests
import bx.intervals.io
#import bx.phylo.newick_tests
#import bx.phylo.phast_tests

suite = unittest.TestSuite( [ #bx.phylo.newick_tests.suite,
                              #bx.phylo.phast_tests.suite,
                              bx.bitset_tests.suite, 
                              bx.align.score_tests.suite,
                              bx.intervals.io.suite ] )
