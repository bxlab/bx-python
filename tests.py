import sys; sys.path.append( "lib" )

import unittest, doctest
import bx.bitset_tests
import bx.align.maf_tests
import bx.align.score_tests
import bx.intervals.io
import bx.phylo.newick_tests
import bx.phylo.phast_tests
import bx.seqmapping_tests
import bx.seq.seq_tests
import bx.seq.nib_tests
import bx.seq.qdna_tests

tests = []

# Modules with doctests
tests.append( doctest.DocTestSuite( bx.intervals.io ) )
tests.append( doctest.DocTestSuite( bx.intervals.intersection ) )

# Standalone test suites
tests.append( bx.phylo.newick_tests.suite )
tests.append( bx.phylo.phast_tests.suite )
tests.append( bx.bitset_tests.suite )
tests.append( bx.align.maf_tests.suite )
tests.append( bx.align.score_tests.suite )
tests.append( bx.seqmapping_tests.suite )
tests.append( bx.seq.seq_tests.suite )
tests.append( bx.seq.nib_tests.suite )
tests.append( bx.seq.qdna_tests.suite )

suite = unittest.TestSuite( tests )
