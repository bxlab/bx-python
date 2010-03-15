import sys, os
import unittest
import tempfile
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
except:
    sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

from bx.arrays.array_tree import ArrayTree, FileArrayTree, FileArrayTreeDict, array_tree_dict_from_reader
from bx.arrays.bed import BedReader
from bx.arrays.wiggle import WiggleReader

class TestArrayTree(unittest.TestCase):
    def setUp(self):
        tree = ArrayTree(10000, 10) # max value of 10000, each block has 10 numbers
        for i in range(5000):
            tree[i] = i
        
        # Insert extra copies to test frequency
        for i in range(3000):
            tree[i] = i
        
        tree.set_range(5000, 9001, 100)
        tree.root.build_summary()
        
        d = {'test': tree}
        f = tempfile.TemporaryFile()
        FileArrayTreeDict.dict_to_file( d, f )
        f.seek(0)
        self.filearraytreedict = FileArrayTreeDict(f)
        self.filearraytree = self.filearraytreedict['test']
        
    def test_get_summary(self):
        f = self.filearraytree
        lvl1 = f.get_summary(0, 1)
        self.assertEqual( map(float, lvl1.sums/lvl1.counts), [4.5, 14.5, 24.5, 34.5, 44.5, 54.5, 64.5, 74.5, 84.5, 94.5])
        lvl2 = f.get_summary(0, 2)
        self.assertEqual( map(float, lvl2.sums/lvl2.counts), [49.5, 149.5, 249.5, 349.5, 449.5, 549.5, 649.5, 749.5, 849.5, 949.5])
        lvl3 = f.get_summary(0, 3)
        self.assertEqual( map(float, lvl3.sums/lvl3.counts), [499.5, 1499.5, 2499.5, 3499.5, 4499.5, 100.0, 100.0, 100.0, 100.0, 100.0])
        lvl2_2 = f.get_summary(3000, 2)
        self.assertEqual( map(float, lvl2_2.sums/lvl2_2.counts), [3049.5, 3149.5, 3249.5, 3349.5, 3449.5, 3549.5, 3649.5, 3749.5, 3849.5, 3949.5])
        
    def test_get_leaf(self):
        f = self.filearraytree
        from_start = [int(i) for i in f.get_leaf(0)]
        from_middle = [int(i) for i in f.get_leaf(5)]
        self.assertEqual(from_start, from_middle)
        self.assertEqual(from_start, range(10))
        
        from_start = [int(i) for i in f.get_leaf(4999)]
        self.assertEqual(from_start, range(4990, 5000))
        
        from_start = [int(i) for i in f.get_leaf(9600)]
        self.assertEqual(from_start, [])
        
    def test_big(self):
        tree = ArrayTree(2147483647, 1000) # What we use for tracks
        for i in range(5000):
            tree[i] = i
        
        # Insert extra copies to test frequency
        for i in range(3000):
            tree[i] = i
        
        tree.set_range(5000, 9001, 100)
        tree.set_range(14000000, 15000000, 200)
        tree.root.build_summary()
        
        d = {'test': tree}
        f = tempfile.TemporaryFile()
        FileArrayTreeDict.dict_to_file( d, f )
        f.seek(0)
        at = FileArrayTreeDict(f)['test']
        
        lvl1 = at.get_summary(14000000, 1)
        avgs = map(float, lvl1.sums/lvl1.counts)
        self.assertEqual( len(avgs), 1000 )
        self.assertEqual( avgs, [ 200 for i in range(0, 1000)] )
    
    
#    def create_bed(self):
#        reader = BedReader( open( "22.bed.txt" ) )
#        temp = tempfile.TemporaryFile()
#        
#        d = array_tree_dict_from_reader( reader, {}, block_size = 1000 )
#
#        for array_tree in d.itervalues():
#            array_tree.root.build_summary()
#
#        FileArrayTreeDict.dict_to_file( d, open("tree.at", "w"), no_leaves=True ) # just summaries
#        
#    def test_bed(self):
#        # self.create_bed()
#        print "bed"
#        at = FileArrayTreeDict( open( "tree.at" ) )['chr22']
#        print map(, at.get_summary(14000000, 1).frequencies)
        
    
    def test_get_frequencies(self):
        f = self.filearraytree
        self.assertEqual( map(float, f.get_summary(0, 1).frequencies), ([20] * 10) )
        self.assertEqual( map(float, f.get_summary(4000, 1).frequencies), ([10] * 10) )
        self.assertEqual( map(float, f.get_summary(0, 2).frequencies), ([200] * 10) )
        self.assertEqual( map(int, f.get_summary(0, 3).frequencies), [2000, 2000, 2000, 1000, 1000, 1000, 1000, 1000, 1000, 1] )
    
    def test_wrong_dictkey(self):
        self.assertRaises(KeyError, self.filearraytreedict.__getitem__, "non-existing")
        
    def test_higher_level_than_tree(self):
        f = self.filearraytree
        self.assertEqual(3, f.levels)
        self.assertRaises(ValueError, f.get_summary, 0, 4)
        

if __name__ == '__main__':
    unittest.main()