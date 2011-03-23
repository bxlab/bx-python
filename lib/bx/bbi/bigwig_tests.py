import sys, os
import unittest
import numpy

try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
except:
    sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

from bx.bbi.bigwig_file import BigWigFile

def allclose( a, b, tol=0.00001 ):
    """
    Like numpy.allclose but treat Nan == Nan 
    """
    d = numpy.absolute( a - b )
    return numpy.all( numpy.isnan( d ) | ( d < tol ) )

class TestBigWig(unittest.TestCase):
    def setUp(self):
        f = open( "test_data/bbi_tests/test.bw" )
        self.bw = BigWigFile(file=f)
        
    def test_get_summary(self):
        data = self.bw.query("chr1", 10000, 20000, 10)
        means = [ x['mean'] for x in data ]
        assert numpy.allclose( map(float, means), [-0.17557571594973645, -0.054009292602539061, -0.056892242431640622, -0.03650328826904297, 0.036112907409667966, 0.0064466032981872557, 0.036949024200439454, 0.076638259887695306, 0.043518108367919923, 0.01554749584197998] )
        
        # Summarize variant
        sd = self.bw.summarize( "chr1", 10000, 20000, 10)
        assert numpy.allclose( sd.sum_data / sd.valid_count, [-0.17557571594973645, -0.054009292602539061, -0.056892242431640622, -0.03650328826904297, 0.036112907409667966, 0.0064466032981872557, 0.036949024200439454, 0.076638259887695306, 0.043518108367919923, 0.01554749584197998] )
        
        # Test min and max for this entire summary region
        data = self.bw.query("chr1", 10000, 20000, 1)
        maxs = [ x['max'] for x in data ]
        mins = [ x['min'] for x in data ]
        self.assertEqual( map(float, maxs), [0.289000004529953] )
        self.assertEqual( map(float, mins), [-3.9100000858306885] )
        
    def test_get_leaf(self):
        data = self.bw.query("chr1", 11000, 11005, 5)
        means = [ x['mean'] for x in data ]
        assert numpy.allclose( map(float, means), [0.050842501223087311, -2.4589500427246094, 0.050842501223087311, 0.050842501223087311, 0.050842501223087311] )
        
        # Test min and max for this entire leaf region
        data = self.bw.query("chr1", 11000, 11005, 1)
        maxs = [ x['max'] for x in data ]
        mins = [ x['min'] for x in data ]
        self.assertEqual( map(float, maxs), [0.050842501223087311] )
        self.assertEqual( map(float, mins), [-2.4589500427246094] )
        
    def test_wrong_nochrom(self):
        data = self.bw.query("chr2", 0, 10000, 10)
        self.assertEqual( data, None )

# Nose test generator
def test_summaries_from_file():
    bw = BigWigFile( file=open( "test_data/bbi_tests/test.bw" ) )
    def check_summary( line ):
        fields = line.split()
        chrom = fields[0]
        start = int( fields[1] )
        end = int( fields[2] )
        n = int( fields[3] )
        t = fields[4]
        values = [ float( v.replace( 'n/a', 'NaN' ) ) for v in fields[5:] ]
        sd = bw.summarize( chrom, start, end, n )
        if t == 'mean':
            print sd.sum_data / sd.valid_count
            print values
            assert allclose( sd.sum_data / sd.valid_count, values )
        elif t == 'min':
            assert allclose( sd.min_val, values )
        elif t == 'max':
            assert allclose( sd.max_val, values )
        #elif t == 'std':
        #    assert numpy.allclose( sd.max_val, values )
    for line in open( "test_data/bbi_tests/test.expectation" ):
        yield check_summary, line

        
if __name__ == '__main__':
    unittest.main()
