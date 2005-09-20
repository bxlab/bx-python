import bx.align.score
import bx.align.maf
import StringIO
import unittest

aligns = [ ( "CCACTAGTTTTTAAATAATCTACTATCAAATAAAAGATTTGTTAATAATAAATTTTAAATCATTAACACTT",
             "CCATTTGGGTTCAAAAATTGATCTATCA----------TGGTGGATTATTATTTAGCCATTAAGGACAAAT", 
             -111 ),
           ( "CCACTAGTTTTTAAATAATCTAC-----AATAAAAGATTTGTTAATAAT---AAATTTTAAATCATTAA-----CACTT",
             "CCATTTGGGTTCAAAAATTGATCTATCA----------TGGTGGAT---TATTATTT-----AGCCATTAAGGACAAAT", 
             -3626 ),
           ( "CCACTAGTTTTTGATTC",
             "CCATTTGGGTTC-----", 
             -299 ),
           ( "CTTAGTTTTTGATCACC",
             "-----CTTGGGTTTACC", 
             -299 ),
           ( "gggaattgaacaatgagaacacatggacacaggaaggggaacatcacacacc----------ggggcctgttgtggggtggggggaag",
             "ggaactagaacaagggagacacatacaaacaacaacaacaacaacacagcccttcccttcaaagagcttatagtctgatggaggagag",
             1690 )
         ]

mafs = """##maf
a score=2883.0
s hg17.chr1             6734 30 + 245522847 CTACCTCAGTGTGGAAGGTGGGCAGTTCTG
s rheMac1.SCAFFOLD71394 9319 30 -     13789 CTACCTCAGTGTGGAAGGTGGGCAGTTCTG

a score=8167.0
s hg17.chr1             41401 40 + 245522847 TGTGTGATTAATGCCTGAGACTGTGTGAAGTAAGAGATGG
s panTro1.chr1          49673 40 + 229575298 TGCGTGATTAATGCCTGAGATTGTGTGAAGTAAAAGATGG
s rheMac1.SCAFFOLD45837 26063 33 -     31516 TGTGTGATTAATGCCTGAGATTGTGTGAAGTAA-------
"""

class BasicTests( unittest.TestCase ):
    def test_scoring_text( self ):
        ss = bx.align.score.hox70
        for t1, t2, score in aligns:
            self.assertEquals( bx.align.score.score_texts( ss, t1, t2 ), score )
            
    def test_align( self ):
        ss = bx.align.score.hox70
        for block in bx.align.maf.Reader( StringIO.StringIO( mafs ) ):
            self.assertEquals( bx.align.score.score_alignment( ss, block ), float( block.score ) )
            
test_classes = [ BasicTests ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )