import unittest

from bx.intervals import *

test_intervals = [ 
    IntervalWithValue( 0,  10,    "1" ),
    IntervalWithValue( 9,   40,   "2" ),
    IntervalWithValue( 5,   11,   "3" ),
    IntervalWithValue( 30,  42,   "4" ),
    IntervalWithValue( 31,  43,   "5" ),
    IntervalWithValue( 11,  1000, "6" ),
    IntervalWithValue( 100, 101,  "7" ),
    IntervalWithValue( 20,  4000, "8" )
]

class TestCase( unittest.TestCase ):

    def testBasic( self ):
        nx = Intersecter()
        for i in test_intervals: nx.add_interval( i )

        r = nx.find( 5, 35 )
        r.sort()

        assert [ '1', '3', '2', '6', '8', '4', '5' ] == [ interval.value for interval in r ]

if __name__ == "__main__":
    unittest.main()
