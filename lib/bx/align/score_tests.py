"""
Tests for `bx.align.score`.
"""

import bx.align.score
import bx.align.maf
import StringIO
import unittest
import sys

from numpy import array, cumsum, allclose

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

nonsymm_scheme = bx.align.score.build_scoring_scheme ( """  A    C    G    T
                                                           91    0  -31 -123
                                                         -114  100 -125  -31
                                                          -31 -125  100 -114
                                                         -123  -31 -114   91 """, 400, 30 )

aligns_for_nonsymm_scheme = [ ( "AAAACCCCGGGGTTTT",
                                "ACGTACGTACGTACGT", 
                                -580 )
                            ]


asymm_scheme = bx.align.score.build_scoring_scheme ( """    01   02    A    C    G    T
                                                       01  200 -200  -50  100  -50  100
                                                       02 -200  200  100  -50  100  -50 """,
                                                       0, 0, gap1='\x00' )

aligns_for_asymm_scheme = [ ( "\x01\x01\x01\x01\x01\x01",
                              "ACGT\x01\x02", 
                              100 )
                          ]


class BasicTests( unittest.TestCase ):

    def test_scoring_text( self ):
        ss = bx.align.score.hox70
        for t1, t2, score in aligns:
            self.assertEquals( bx.align.score.score_texts( ss, t1, t2 ), score )
            
    def test_align( self ):
        ss = bx.align.score.hox70
        for block in bx.align.maf.Reader( StringIO.StringIO( mafs ) ):
            self.assertEquals( bx.align.score.score_alignment( ss, block ), float( block.score ) )
            
    def test_accumulate( self ):
        ss = bx.align.score.hox70
        self.assert_( allclose( bx.align.score.accumulate_scores( ss, "-----CTTT", "CTTAGTTTA"  ),
                           cumsum( array( [ -430, -30, -30, -30, -30, -31, 91, 91, -123 ] ) ) ) )
        self.assert_( allclose( bx.align.score.accumulate_scores( ss, "-----CTTT", "CTTAGTTTA", skip_ref_gaps=True ),
                           cumsum( array( [ -581, 91, 91, -123 ] ) ) ) )

    def test_nonsymm_scoring( self ):
        ss = nonsymm_scheme
        for t1, t2, score in aligns_for_nonsymm_scheme:
            self.assertEquals( bx.align.score.score_texts( ss, t1, t2 ), score )

    def test_asymm_scoring( self ):
        ss = asymm_scheme
        for t1, t2, score in aligns_for_asymm_scheme:
            self.assertEquals( bx.align.score.score_texts( ss, t1, t2 ), score )
   
test_classes = [ BasicTests ]
suite = unittest.TestSuite( [ unittest.makeSuite( c ) for c in test_classes ] )
