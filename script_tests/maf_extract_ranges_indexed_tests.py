import unittest

import base


class Test(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_extract_ranges_indexed.py ./test_data/maf_tests/mm8_chr7_tiny.maf -c -m 5 -p mm8."
    input_stdin = base.TestFile(filename="./test_data/maf_tests/dcking_ghp074.bed")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/dcking_ghp074.maf")


class TestAccessNotRef(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_extract_ranges_indexed.py ./test_data/maf_tests/mm8_chr7_tiny.maf -c -m 5 -p hg18."
    input_stdin = base.TestFile(filename="./test_data/maf_tests/hg18.bed")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/test_hg18.maf")


class TestAccessRef(base.BaseScriptTest, unittest.TestCase):
    command_line = (
        "./scripts/maf_extract_ranges_indexed.py ./test_data/maf_tests/mm8_chr7_tiny_mm8_ind.maf -c -m 5 -p mm8."
    )
    input_stdin = base.TestFile(filename="./test_data/maf_tests/dcking_ghp074.bed")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/dcking_ghp074.maf")


class TestAccessNotRefNotIndexed(base.BaseScriptTest, unittest.TestCase):
    command_line = (
        "./scripts/maf_extract_ranges_indexed.py ./test_data/maf_tests/mm8_chr7_tiny_mm8_ind.maf -c -m 5 -p hg18."
    )
    input_stdin = base.TestFile(filename="./test_data/maf_tests/hg18.bed")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/empty.maf")


class TestELines(base.BaseScriptTest, unittest.TestCase):
    command_line = (
        "./scripts/maf_extract_ranges_indexed.py ./test_data/maf_tests/mm10_chr12_lessspe.maf -c -m 5 -p mm10."
    )
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm10_chr12.bed")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_slice.maf")
