import sys, os
import unittest
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
except:
    sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

# from bx.intervals.cluster import ClusterTree
from cluster import ClusterTree

class TestCluster(unittest.TestCase):
    def setUp(self):
        self.tree = ClusterTree(0, 0)
    
    def insertpairs(self, pairs):
        for i, (s, e) in enumerate(pairs):
            self.tree.insert(s, e, i)
    
    def test_merge_case(self):
        pairs = [(3, 4), (6, 7), (9, 10), (1, 2), (3, 8)]
        self.insertpairs(pairs)
        
        self.assertEqual( [(1, 2, [3]), (3, 8, [0, 1, 4]), (9, 10, [2])], self.tree.getregions() )
    
    def test_trivial(self):
        pairs = [ (1, 4), (4, 5) ]
        self.insertpairs(pairs)
        
        self.assertEqual( [(1, 5, [0, 1])], self.tree.getregions() )
        
    def test_easymerge(self):
        pairs = [ (1, 2), (4, 5), (2, 4) ]
        self.insertpairs(pairs)

        self.assertEqual( [(1, 5, [0, 1, 2])], self.tree.getregions() )
        
    def test_hardmerge(self):
        pairs = [ (1, 2), (8, 9), (3, 4), (5, 6), (7, 8), (1, 10) ]
        self.insertpairs(pairs)

        self.assertEqual( [(1, 10, [0, 1, 2, 3, 4, 5])], self.tree.getregions() )
    
    def test_duplicates(self):
        pairs = [ (1, 1), (1, 2), (3, 4), (3, 4), (1, 4) ]
        self.insertpairs(pairs)

        self.assertEqual( [(1, 4, [0, 1, 2, 3, 4])], self.tree.getregions() )
        
    def test_startbeforeend(self):
        self.assertRaises(ValueError, self.tree.insert, 4, 2, 0)
        
    def test_large_sorted(self):
        upto = 100000
        pairs = [ (2*i + 1, 2*i + 2) for i in range(upto) ]
        self.insertpairs(pairs)
        self.tree.insert( 0, upto*3, upto )
        self.assertEqual( [ (0, upto*3, [x for x in range(upto+1)]) ], self.tree.getregions() )
        
    def test_minregions(self):
        self.tree = ClusterTree(0, 2)
        pairs = [(3, 4), (6, 7), (9, 10), (1, 2), (3, 8)]
        self.insertpairs(pairs)
        
        self.assertEqual( [(3, 8, [0, 1, 4])], self.tree.getregions() )
        
    def test_distance(self):
        self.tree = ClusterTree(1, 0)
        pairs = [(3, 4), (6, 7), (9, 10), (1, 2), (3, 8)]
        self.insertpairs(pairs)
        
        self.assertEqual( [(1, 10, [0, 1, 2, 3, 4])], self.tree.getregions() )
        
    def test_merge_left_right(self):
        pairs = [(6, 7, 1), (1, 2, 3), (9, 10, 2), (3, 4, 0), (3, 8, 4)]
        for s, e, i in pairs:
            self.tree.insert(s, e, i)

        self.assertEqual( [(1, 2, [3]), (3, 8, [0, 1, 4]), (9, 10, [2])], self.tree.getregions() )
    
    def test_larger(self):
        pairs = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12), (13, 14), (15, 16), (17, 18), (19, 20),
                (1, 3), (4, 10), (10, 15), (15, 20), (21, 22)]
        self.insertpairs(pairs)
        
        self.assertEqual( [(1, 20, [x for x in range(14)]), (21, 22, [14])], self.tree.getregions() )
    
    def test_another(self):
        pairs = [(3, 4, 1), (13, 14, 6), (21, 22, 14), (5, 6, 2), (4, 10, 11), (1, 2, 0), (11, 12, 5), (1, 3, 10), (7, 8, 3), (15, 16, 7), (15, 20, 13), (19, 20, 9), (10, 15, 12), (17, 18, 8), (9, 10, 4)]
        # pairs = [(3, 4, 1), (13, 14, 6), (21, 22, 14), (5, 6, 2), (4, 10, 11), (1, 2, 0), (11, 12, 5), (1, 3, 10), (7, 8, 3), (15, 16, 7), (15, 20, 13), (19, 20, 9), (10, 15, 12), (9, 10, 4)]
        for s, e, i in pairs:
            self.tree.insert(s, e, i)
        
        self.assertEqual( [(1, 20, [x for x in range(14)]), (21, 22, [14])], self.tree.getregions() )
        
    def test_none(self):
        pairs = []
        self.insertpairs(pairs)

        self.assertEqual( [], self.tree.getregions() )
        

if __name__ == '__main__':
    unittest.main()