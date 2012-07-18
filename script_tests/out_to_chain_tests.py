import base
import unittest

class Test( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/out_to_chain.py ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.out --chrsizes ./test_data/epo_tests/hg19.chrom.sizes ./test_data/epo_tests/mm9.chrom.sizes"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain" )

unittest.main()
