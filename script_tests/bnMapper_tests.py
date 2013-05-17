import base
import unittest

class Test1( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/out_to_chain.py ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.out --chrsizes ./test_data/epo_tests/hg19.chrom.sizes ./test_data/epo_tests/mm9.chrom.sizes"
    command_line = "./scripts/bnMapper.py ./test_data/epo_tests/hpeaks.bed ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hpeaks.mapped.bed4" )

class Test2( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/bnMapper.py -fBED12 ./test_data/epo_tests/hpeaks.bed ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hpeaks.mapped.bed12" )

class Test3( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/bnMapper.py -g9 ./test_data/epo_tests/hpeaks.bed ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hpeaks.mapped.bed4" )

class Test4( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/bnMapper.py -g3 ./test_data/epo_tests/hpeaks.bed ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hpeaks.mapped.nopeak2.bed4" )

class Test5( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/bnMapper.py -g9 -t0.67 ./test_data/epo_tests/hpeaks.bed ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hpeaks.mapped.bed4" )

class Test6( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/bnMapper.py -g9 -t0.7 ./test_data/epo_tests/hpeaks.bed ./test_data/epo_tests/epo_547_hs_mm_12way_mammals_65.chain"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hpeaks.mapped.nopeak2.bed4" )

class Test6( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/bnMapper.py ./test_data/epo_tests/hg19_one_peak.bed ./test_data/epo_tests/hg19.mm9.rBest.chain.gz"
    output_stdout = base.TestFile( filename="./test_data/epo_tests/hg19_one_peak.mapped.bed" )

unittest.main()

