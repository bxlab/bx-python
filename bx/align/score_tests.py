import bx.align.score
import unittest

aligns = [ ( "CCACTAGTTTTTAAATAATCTACTATCAAATAAAAGATTTGTTAATAATAAATTTTAAATCATTAACACTT",
             "CCATTTGGGTTCAAAAATTGATCTATCA----------TGGTGGATTATTATTTAGCCATTAAGGACAAAT", 
             -111 ),
           ( "gggaattgaacaatgagaacacatggacacaggaaggggaacatcacacacc----------ggggcctgttgtggggtggggggaag",
             "ggaactagaacaagggagacacatacaaacaacaacaacaacaacacagcccttcccttcaaagagcttatagtctgatggaggagag",
             1690 )   
         ]

class BasicTests( unittest.TestCase ):
    def test_scoring( self ):
        ss = bx.align.score.hox70
        for t1, t2, score in aligns:
            self.assertEquals( bx.align.score.score( ss, t1, t2 ), score )
            
test_classes = [ BasicTests ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )