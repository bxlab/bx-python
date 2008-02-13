import base
import unittest

class Test( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/maf_extract_ranges_indexed.py /depot/data2/galaxy/mm8/align/multiz17way/*.maf.lzo -c -p mm8. < ${bed}" 
    input_bed = base.TestFile( filename="./test_data/maf_tests/dcking_ghp074.bed" )
    output_stdout = base.TestFile( filename="./test_data/maf_tests/dcking_ghp074.maf" )
