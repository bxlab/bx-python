import unittest

import base


class Test(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_select.py ${features}"
    input_features = base.TestFile("""0
                                       0
                                       0
                                       0
                                       0
                                       0
                                       0
                                       1""")
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_no_index.maf")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_last_selected.maf")


class TestWithE(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_select.py ${features}"
    input_features = base.TestFile("""0
                                       1
                                       0
                                       0
                                       0
                                       0
                                       0
                                       0""")
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe.maf")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe_one_selected.maf")
