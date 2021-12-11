import unittest

import base


class Test1(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_build_index.py ${maf} ${maf_index}"
    input_maf = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf")
    output_maf_index = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.index")


class Test2(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_build_index.py ${maf_bz2} ${maf_index}"
    input_maf_bz2 = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.bz2")
    output_maf_index = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.index")


class Test3(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_build_index.py ${maf_lzo} ${maf_index}"
    input_maf_lzo = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.lzo")
    output_maf_index = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.index")


class TestindexOnlyRef(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_build_index.py -s mm8 ${maf} ${maf_index}"
    input_maf = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_mm8_ind.maf")
    output_maf_index = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_mm8_ind.maf.index")


class TestindexWithElines(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_build_index.py ${maf} ${maf_index}"
    input_maf = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe.maf")
    output_maf_index = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe.maf.index")
