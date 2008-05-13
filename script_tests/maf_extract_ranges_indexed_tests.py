import base
import unittest

class Test( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/maf_extract_ranges_indexed.py ./test_data/maf_tests/mm8_chr7_tiny.maf -c -p mm8. < ${bed}" 
    input_bed = base.TestFile( filename="./test_data/maf_tests/dcking_ghp074.bed" )
    output_stdout = base.TestFile( filename="./test_data/maf_tests/dcking_ghp074.maf" )
