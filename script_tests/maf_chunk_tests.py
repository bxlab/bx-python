import unittest

import base


class Test1(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_chunk.py 1 ${out_dir}"
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf")
    out_dir = "./test_data/maf_tests/chunk1"


class Test2(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_chunk.py 1000 ${out_dir}"
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf")
    out_dir = "./test_data/maf_tests/chunk1000"
