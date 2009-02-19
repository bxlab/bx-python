import sys, os
import unittest
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
except:
    sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalNode

class NeighborTestCase(unittest.TestCase):

    def setUp(self):
        iv = IntervalNode(Interval(50, 59))
        for i in range(0, 110, 10):
            if i == 50: continue
            f = Interval(i, i + 9)
            iv = iv.insert(f)
        self.intervals = iv

    def test_left(self):
        iv = self.intervals 
        self.assertEqual(str(iv.left(Interval(60, 70), n=2)), str([Interval(50, 59), Interval(40, 49)]))

        for i in range(10, 100, 10):
            f = Interval(i, i)
            r = iv.left(f, max_dist=10, n=1)
            self.assertEqual(r[0].end,  i - 1)

    def test_toomany(self):
        iv = self.intervals 
        self.assertEqual(len(iv.left(Interval(60, 70), n=200)) , 6)


    def test_right(self):
        iv = self.intervals 
        self.assertEqual(str(iv.left(Interval(60, 70), n=2)), str([Interval(50, 59), Interval(40, 49)]))

        def get_right_start(b10):
            r = iv.right(Interval(b10, b10 + 1), n=1)
            assert len(r) == 1
            return r[0].start
        
        for i in range(10, 100, 10):
            self.assertEqual(get_right_start(i), i + 10)

        for i in range(0, 100, 10):
            f = Interval(i - 1, i - 1)
            r = iv.right(f, max_dist=10, n=1)
            self.assertEqual(r[0].start, i)

    def test_upstream(self):
        iv = self.intervals 
        upstreams = iv.upstream(Interval(59, 60), n=200)
        for u in upstreams:
            self.assertTrue(u.end < 59)

        upstreams = iv.upstream(Interval(60, 70, strand=-1), n=200)
        for u in upstreams:
            self.assertTrue(u.start > 70)


        upstreams = iv.upstream(Interval(58, 58, strand=-1), n=200)
        for u in upstreams:
            self.assertTrue(u.start > 59)


    def test_downstream(self):
        iv = self.intervals 
        downstreams = iv.downstream(Interval(59, 60), n=200)
        for d in downstreams:
            self.assertTrue(d.start > 60)

        downstreams = iv.downstream(Interval(59, 60, strand=-1), n=200)
        for d in downstreams:
            self.assertTrue(d.start < 59)
        

    def test_n(self):
        iv = self.intervals
        for i in range(0, 90, 10):
            f = Interval(i + 1, i + 1)
            r = iv.right(f, max_dist=20, n=2)
            self.assertEqual(r[0].start, i + 10)
            self.assertEqual(r[1].start, i + 20)


class LotsaTestCase(unittest.TestCase):
    """ put lotsa data in the tree and make sure it works"""
    def setUp(self):
        iv = IntervalNode(Interval(1, 2))
        self.max = 1000000
        for i in range(0, self.max, 10):
            f = Interval(i, i)
            iv = iv.insert(f)

        for i in range(600):
            iv = iv.insert(Interval(0, 1))
        self.intervals = iv



    def test_count(self):
        iv = self.intervals
        
        r = iv.right(Interval(1, 1), n=33)
        self.assertEqual(len(r), 33)

        l = iv.left(Interval(1, 1), n=33)
        self.assertEqual(len(l), 1)

        u = iv.upstream(Interval(1, 1, strand=-1), n=9999)
        self.assertEqual(len(u), 250)

        # now increase max_dist
        u = iv.upstream(Interval(1, 1, strand=-1), n=9999, max_dist=99999)
        self.assertEqual(len(u), 9999)


    def test_max_dist(self):        
        iv = self.intervals
        r = iv.right(Interval(1, 1), max_dist=0, n=10)
        self.assertEqual(len(r), 0)

        for n, d in enumerate(range(10, 1000, 10)):
            r = iv.right(Interval(1, 1), max_dist=d, n=10000)
            self.assertEqual(len(r), n + 1) 

    def test_find(self):
        iv = self.intervals
        random =  __import__("random")
        for t in range(25):
            start = random.randint(0, self.max - 10000)
            end  = start + random.randint(100, 10000)

            results = iv.find(start, end)
            for feat in results:
                self.assertTrue(
                        (feat.end >= start and feat.end <= end) 
                            or 
                        (feat.start <= end and feat.start >= start)
                        )


if __name__ == "__main__":
    unittest.main()
