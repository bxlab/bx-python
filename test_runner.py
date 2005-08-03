#!/usr/bin/env python

import unittest
import sys
import distutils.util

lib_dir = "./build/lib.%s-%s" % ( distutils.util.get_platform(), sys.version[0:3] ) 
sys.path = [ lib_dir ] + sys.path

import bx.bitset_tests

tests = [ bx.bitset_tests.suite ]

if __name__ == "__main__": 
    test_runner = unittest.TextTestRunner( verbosity=2 )
    result = test_runner.run( unittest.TestSuite( tests ) )
    sys.exit(not result.wasSuccessful())
