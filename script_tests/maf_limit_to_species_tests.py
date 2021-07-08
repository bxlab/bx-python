import unittest

import base


class Test1(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_limit_to_species.py mm8,rn4"
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_no_index.maf")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_only_mouse_rat.maf")


class TestWithE(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_limit_to_species.py mm10,bosTau7,loxAfr3"
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe.maf")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe_onlymouse_cow_elephant.maf")
