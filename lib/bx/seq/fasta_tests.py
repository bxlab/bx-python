import unittest
import os.path
import fasta

# Same sequence data as stored in test.fa

test_seq = "TGGAGGCATTTGTGATTCAATAGATGCAGAAAGAAACCTTCCTAGAGCTG" \
         + "GCGTTCTCTAACTAAAAGTGGAAAGTTCTGAGGAATGAGGACTGTTATAA" \
         + "ATCCCACCCCACACCGCACCTTCTCCAGGGAAGTTTCATGGCCGTGAAGA" \
         + "GGACAGAAAGTGAGAACCAAGATggaactgaataaacaagcttcacactg" \
         + "ttagtttccccatatgcttaccttcccacagatgccaaccttggaggcct" \
         + "aagaggcctagaatattatcctttgtctgatcatttctctacaaatttat" \
         + "tgttctttgttaagatgctacataagcccaaattctaaccacccctttga" \
         + "gttacccatcatcaagtttctcccatgtg"

test_seq_len = len( test_seq )

class FASTATestCase( unittest.TestCase ):

    def test_get( self ):
        fastafile = fasta.FastaFile( file( os.path.join('lib','bx','seq','test.fa'), "rb" ) )
        do_test_get(fastafile, 0, test_seq_len)
        do_test_get(fastafile, 0, 40)
        do_test_get(fastafile, test_seq_len - 40, 40)

def do_test_get( fastafile, start, len ):
    assert fastafile.get( start, len ) == test_seq[start:start+len]

test_classes = [ FASTATestCase ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
