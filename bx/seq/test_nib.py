import unittest
import nib

# Same sequence data as stored in test.nib

test_seq = "TGGAGGCATTTGTGATTCAATAGATGCAGAAAGAAACCTTCCTAGAGCTGGCGTTCTCTAACTAAAAGTGGAAAGTTCTGAGGAATGAGGACTGTTATAAATCCCACCCCACACCGCACCTTCTCCAGGGAAGTTTCATGGCCGTGAAGAGGACAGAAAGTGAGAACCAAGATggaactgaataaacaagcttcacactgttagtttccccatatgcttaccttcccacagatgccaaccttggaggcctaagaggcctagaatattatcctttgtctgatcatttctctacaaatttattgttctttgttaagatgctacataagcccaaattctaaccacccctttgagttacccatcatcaagtttctcccatgtg"

test_seq_len = len( test_seq )

class NIBTestCase( unittest.TestCase ):

    def test_get( self ):
        nibfile = nib.NibFile( file( "test.nib" ) )
        # Try all combinations of even / odd boundaries
        do_test_get( nibfile, 0, 10 )
        do_test_get( nibfile, 1, 10 )
        do_test_get( nibfile, 0, 11 )
        do_test_get( nibfile, 1, 11 )
        # Test near end of file also
        do_test_get( nibfile, test_seq_len - 10, 10 )
        do_test_get( nibfile, test_seq_len - 11, 11 )

def do_test_get( nibfile, start, len ):
    assert nibfile.get( start, len ) == test_seq[start:start+len]

if __name__ == "__main__": unittest.main()
