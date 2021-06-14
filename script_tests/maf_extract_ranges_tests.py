import unittest

import base


class Test(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_extract_ranges.py ${inverval_file} 0"
    input_inverval_file = base.TestFile("80082367 80083066")
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_no_index.maf")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/dcking_ghp074.maf")

class TestElines(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_extract_ranges.py ${inverval_file} 0 -m 5"
    input_inverval_file = base.TestFile("56694985 56695040")
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe.maf")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_slice2.maf")
