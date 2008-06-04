import pwm
from numpy import allclose, isnan

def test_create():
    m = pwm.FrequencyMatrix.from_rows( ['A','C','G','T'], get_ctcf_rows() )
    # Alphabet sort
    assert m.sorted_alphabet == [ 'A', 'C', 'G', 'T' ]
    # Character to index mapping
    assert m.char_to_index[ ord('A') ] == 0
    assert m.char_to_index[ ord('C') ] == 1
    assert m.char_to_index[ ord('G') ] == 2
    assert m.char_to_index[ ord('T') ] == 3
    assert m.char_to_index[ ord('Q') ] == -1
    # Values
    assert allclose( m.values[0],  [ 2620, 2052, 3013, 2314 ] )
    assert allclose( m.values[19],  [ 3144, 3231, 3056, 567 ] )
    
def test_scoring():
    m = pwm.FrequencyMatrix.from_rows( ['A','C','G','T'], get_ctcf_rows() )
    # Stormo method
    sm = m.to_stormo_scoring_matrix()
    # Forward matches
    assert allclose( sm.score_string( "AATCACCACCTCCTGGCAGG" )[0], -156.8261261 )
    assert allclose( sm.score_string( "TGCCTGCCTCTGTAGGCTCC" )[0], -128.8106842 )
    assert allclose( sm.score_string( "GTTGCCAGTTGGGGGAAGCA" )[0], 4.65049839 )
    assert allclose( sm.score_string( "GCAGACACCAGGTGGTTCAG" )[0], 1.60168743 )
    # Reverse matches
    rc = sm.reverse_complement()
    assert allclose( rc.score_string( "AATCACCACCTCCTGGCAGG" )[0], 0.014178276062 )
    assert allclose( rc.score_string( "TGCCTGCCTCTGTAGGCTCC" )[0], 0.723828315735 )
    assert allclose( rc.score_string( "GTTGCCAGTTGGGGGAAGCA" )[0], -126.99407196 )
    assert allclose( rc.score_string( "GCAGACACCAGGTGGTTCAG" )[0], -86.9560623169 )
    # Nothing valid
    assert isnan( sm.score_string_with_gaps( "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" ) ).all()
    # Too short
    assert isnan( sm.score_string( "TTTT" ) ).all()

def test_scoring_with_gaps():
    m = pwm.FrequencyMatrix.from_rows( ['A','C','G','T'], get_ctcf_rows() )
    # Stormo method
    sm = m.to_stormo_scoring_matrix()
    # Forward matches
    assert allclose( sm.score_string_with_gaps( "GTTGCCAGT----TGGGGGAAGCATTT---AA" )[0], 4.65049839 )
    assert allclose( sm.score_string_with_gaps( "GCAGA--CACCAGGTGG--TTCAG---" )[0], 1.60168743 )
    assert allclose( sm.score_string_with_gaps( "----GTTGCCAGTTGGGGGAAGCA" )[4], 4.65049839 )
    assert allclose( sm.score_string_with_gaps( "TTT--GTT--GCCA--GTTGGGG-G-A-A-G-C-A-" )[5], 4.65049839 )
    assert isnan( sm.score_string_with_gaps( "TTT--GTT--GCCA--GTTGGGG-G-A-A-G-C-A-" )[4] )
    # Nothing valid
    assert isnan( sm.score_string_with_gaps( "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" ) ).all()
    assert isnan( sm.score_string_with_gaps( "------------------------------------" ) ).all()
    # Too short
    assert isnan( sm.score_string_with_gaps( "TTTT" ) ).all()
    assert isnan( sm.score_string_with_gaps( "TTTT----" ) ).all()


def get_ctcf_rows():
    """
    The CTCF primary site motif
    """
    return [
         [   2620 ,  2052  ,  3013  ,  2314   ],
         [      0 ,  3580  ,  1746  ,  4672   ],
         [   2008 ,  1790  ,  4497  ,  1703   ],
         [   3362 ,     0  ,  6637  ,     0   ],
         [      0 , 10000  ,     0  ,     0   ],
         [      0 , 10000  ,     0  ,     0   ],
         [   7467 ,     0  ,  1310  ,  1222   ],
         [    786 ,  4890  ,  4323  ,     0   ],
         [   1179 ,  6288  ,   829  ,  1703   ],
         [  10000 ,     0  ,     0  ,     0   ],
         [      0 ,     0  , 10000  ,     0   ],
         [   4847 ,     0  ,  5152  ,     0   ],
         [      0 ,     0  ,  6200  ,  3799   ],
         [      0 ,     0  , 10000  ,     0   ],
         [      0 ,     0  , 10000  ,     0   ],
         [   1572 ,  7467  ,     0  ,   960   ],
         [   3842 ,     0  ,  5545  ,   611   ],
         [      0 ,  5895  ,  4104  ,     0   ],
         [   1615 ,  4192  ,  1397  ,  2794   ],
         [   3144 ,  3231  ,  3056  ,   567   ]
    ]