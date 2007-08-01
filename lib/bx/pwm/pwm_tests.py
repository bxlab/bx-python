import unittest
import sys
import bx.pwm.position_weight_matrix as pwm
from StringIO import StringIO

basicPwm = \
""">MA0101 c-REL REL
0   5   8   4
0   1   15  1
1   0   15  1
5   1   9   2
6   5   3   3
5   1   1   10
1   0   0   16
2   0   0   15
0   15  0   2
1   16  0   0
"""

transfacPwm = \
"""ID  TATA
XX
P0    A    C    G    T
01   33   73   78   16      S
02   10   24   11  155      T
03  176    3    2   19      A
04    2    7    3  188      T
05  178    2    3   17      A
06  133    2    2   63      A
07  183    3   10    4      A
08  112    2   24   62      W
09   78   26   80   16      R
10   29   72   75   24      N
11   42   74   68   16      N
12   42   65   66   27      N
13   41   60   67   32      N
14   35   54   72   39      N
15   40   51   73   36      N
XX
"""

background = { 'A':.28,'C':.21, 'G':.24, 'T':.27 }

dSeq = "ACCGAGTTAGCGTAAA"
dScoresExpected = "-15.3697 0.4240 -16.5309 0.4027"

qSeq = [{'A':0.27,'C':0.34,'G':0.07,'T':0.32},
        {'A':0.24,'C':0.32,'G':0.09,'T':0.35},
        {'A':0.80,'C':0.11,'G':0.03,'T':0.06},
        {'A':0.07,'C':0.22,'G':0.37,'T':0.34},
        {'A':0.07,'C':0.44,'G':0.03,'T':0.46},
        {'A':0.43,'C':0.04,'G':0.18,'T':0.35},
        {'A':0.84,'C':0.14,'G':0.01,'T':0.01},
        {'A':0.31,'C':0.52,'G':0.13,'T':0.04},
        {'A':0.22,'C':0.22,'G':0.45,'T':0.11},
        {'A':0.36,'C':0.15,'G':0.42,'T':0.07},
        {'A':0.11,'C':0.78,'G':0.07,'T':0.04},
        {'A':0.07,'C':0.16,'G':0.64,'T':0.13},
        {'A':0.34,'C':0.59,'G':0.03,'T':0.04},
        {'A':0.32,'C':0.15,'G':0.07,'T':0.46},
        {'A':0.07,'C':0.03,'G':0.59,'T':0.31}]

qScoresExpected = "4.1106 0.7810"

class PWMTestCase (unittest.TestCase):

    def testReader(self):

        # test basic format: i.e. for jaspar
        wms = [wm for wm in pwm.Reader(StringIO(basicPwm),format="basic", \
                          background=background,score_correction=False)]
        assert len(wms) == 1

        # test transfac format
        wms = [wm for wm in pwm.Reader(StringIO(transfacPwm),format="transfac", \
                          background=background,score_correction=False)]
        assert len(wms) == 1

        wm = wms[0]
        dScores = wm.score_seq(dSeq)
        assert len(dScores) == 2
        assert "%.4f %.4f %.4f %.4f" % (dScores[0][0],dScores[0][1],dScores[1][0],dScores[1][1]) == dScoresExpected

        qdSeq = []
        for (ix,nt) in enumerate(dSeq):
            qdSeq.append(dict())
            qdSeq[ix][nt] = 1.0
        qScores = wm.score_seq(qdSeq)
        assert len(qScores) == 2
        assert "%.4f %.4f %.4f %.4f" % (qScores[0][0],qScores[0][1],qScores[1][0],qScores[1][1]) == dScoresExpected

        qScores = wm.score_seq(qSeq)
        assert len(qScores) == 1
        assert "%.4f %.4f" % (qScores[0][0],qScores[0][1]) == qScoresExpected

test_classes = [PWMTestCase]
suite = unittest.TestSuite ([unittest.makeSuite (c) for c in test_classes])
