from Numeric import *

class ScoringScheme( object ):
    def __init__( self, gap_open, gap_extend, input_size=128, default=-100, typecode=Int ):
        self.table = zeros( (input_size, input_size), typecode=typecode )
        self.table *= default
        self.gap_open = gap_open
        self.gap_extend = gap_extend
    def set_score( self, a, b, val ):
        self.table[a,b] = val
    
def build_scoring_scheme( s, gap_open, gap_extend ):
    """
    Initialize scoring scheme from a blastz style text blob, first line
    specifies the bases for each row/col, subsequent lines contain the
    corresponding scores.
    """
    ss = ScoringScheme( gap_open, gap_extend )
    lines = s.split( "\n" )
    chars = lines[0].split()
    for i, line in enumerate( lines[1:] ):
        for j, score in enumerate( map( int, line.split() ) ):
            ss.set_score( ord( chars[i].lower() ), ord( chars[j].lower() ), score )
            ss.set_score( ord( chars[i].upper() ), ord( chars[j].upper() ), score )
    return ss
            
def score_alignment( scoring_scheme, a ):
    score = 0
    ncomps = len( a.components )
    for i in range( ncomps ):
        for j in range( i+1, ncomps ):
            score += score_texts( scoring_scheme, a.components[i].text, a.components[j].text )
    return score
    

def score_texts( scoring_scheme, text1, text2 ):
    rval = 0
    last_gap = False
    for i in range( len( text1 ) ):
        a = text1[i]
        b = text2[i]
        if a == '-' or b == '-':
            rval -= scoring_scheme.gap_extend
            if not last_gap:
                rval -= scoring_scheme.gap_open
                last_gap = True
        else:   
            rval += scoring_scheme.table[ord(a),ord(b)]
            last_gap = False
    return rval
            
hox70 = build_scoring_scheme( """  A    C    G    T
                                  91 -114  -31 -123
                                -114  100 -125  -31
                                 -31 -125  100 -114
                                -123  -31 -114   91 """, 400, 30 )