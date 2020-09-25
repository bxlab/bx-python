import unittest

import base


class Test(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/lzop_build_offset_table.py"
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.lzo")
    output_stdout = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny.maf.lzot")
