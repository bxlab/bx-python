"""
Tests for `bx.phylo.newick`.
"""

from bx.phylo.newick import *
import unittest

trees = [ r"(B:6.0,(A:5.0,C:3.0,'Foo ''bar':4.0)Q_X:5.0,D:11.0)label;",
          "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,(( monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);",
          "(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460);",
          "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);",
          "(B,(A,C,E),D);",
          "(,(,,),);",
          "(A,(B,C),D);",
          "((A,D),(C,B));"
]

results = [ ( Tree( 'label', [Edge( 6.0, Tree( 'B', None ) ), Edge( 5.0, Tree( 'Q X', [Edge( 5.0, Tree( 'A', None ) ), Edge( 3.0, Tree( 'C', None ) ), Edge( 4.0, Tree( "Foo 'bar", None ) )] ) ), Edge( 11.0, Tree( 'D', None ) )] ) ),
            ( Tree( None, [Edge( 0.84599999999999997, Tree( None, [Edge( 19.199590000000001, Tree( 'raccoon', None ) ), Edge( 6.8004100000000003, Tree( 'bear', None ) )] ) ), Edge( 3.8738199999999998, Tree( None, [Edge( 7.5297299999999998, Tree( None, [Edge( 11.997, Tree( 'sea lion', None ) ), Edge( 12.003, Tree( 'seal', None ) )] ) ), Edge( 2.0945999999999998, Tree( None, [Edge( 20.592009999999998, Tree( None, [Edge( 100.8593, Tree( 'monkey', None ) ), Edge( 47.140689999999999, Tree( 'cat', None ) )] ) ), Edge( 18.879529999999999, Tree( 'weasel', None ) )] ) )] ) ), Edge( 25.461539999999999, Tree( 'dog', None ) )] ) ),
            ( Tree( None, [Edge( 0.69394999999999996, Tree( 'Bovine', None ) ), Edge( 0.54939000000000004, Tree( None, [Edge( 0.36079, Tree( 'Gibbon', None ) ), Edge( 0.15057000000000001, Tree( None, [Edge( 0.33635999999999999, Tree( 'Orang', None ) ), Edge( 0.061240000000000003, Tree( None, [Edge( 0.17147000000000001, Tree( 'Gorilla', None ) ), Edge( 0.083860000000000004, Tree( None, [Edge( 0.19267999999999999, Tree( 'Chimp', None ) ), Edge( 0.11927, Tree( 'Human', None ) )] ) )] ) )] ) )] ) ), Edge( 1.2145999999999999, Tree( 'Mouse', None ) )] ) ),
            ( Tree( None, [Edge( 0.69394999999999996, Tree( 'Bovine', None ) ), Edge( 0.54939000000000004, Tree( None, [Edge( 0.36079, Tree( 'Hylobates', None ) ), Edge( 0.15057000000000001, Tree( None, [Edge( 0.33635999999999999, Tree( 'Pongo', None ) ), Edge( 0.061240000000000003, Tree( None, [Edge( 0.17147000000000001, Tree( 'G. Gorilla', None ) ), Edge( 0.083860000000000004, Tree( None, [Edge( 0.19267999999999999, Tree( 'P. paniscus', None ) ), Edge( 0.11927, Tree( 'H. sapiens', None ) )] ) )] ) )] ) )] ) ), Edge( 1.2145999999999999, Tree( 'Rodent', None ) )] ) ),
            ( Tree( None, [Edge( None, Tree( 'B', None ) ), Edge( None, Tree( None, [Edge( None, Tree( 'A', None ) ), Edge( None, Tree( 'C', None ) ), Edge( None, Tree( 'E', None ) )] ) ), Edge( None, Tree( 'D', None ) )] ) ),
            ( Tree( None, [Edge( None, Tree( None, None ) ), Edge( None, Tree( None, [Edge( None, Tree( None, None ) ), Edge( None, Tree( None, None ) ), Edge( None, Tree( None, None ) )] ) ), Edge( None, Tree( None, None ) )] ) ),
            ( Tree( None, [Edge( None, Tree( 'A', None ) ), Edge( None, Tree( None, [Edge( None, Tree( 'B', None ) ), Edge( None, Tree( 'C', None ) )] ) ), Edge( None, Tree( 'D', None ) )] ) ),
            ( Tree( None, [Edge( None, Tree( None, [Edge( None, Tree( 'A', None ) ), Edge( None, Tree( 'D', None ) )] ) ), Edge( None, Tree( None, [Edge( None, Tree( 'C', None ) ), Edge( None, Tree( 'B', None ) )] ) )] ) ),
            ]

def tests(): 
    for i in range(len(trees)):
        def _( s, r ):
            assert newick_parser.parse_string( s ) == r
        _.description = "check tree parsing " + str(i)
        yield _, trees[i], results[i] 
