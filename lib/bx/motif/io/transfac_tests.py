from StringIO import StringIO
import transfac
from numpy import allclose

sample = """
VV  TRANSFAC MATRIX TABLE, Rel.3.2 26-06-1997
XX
//
AC  a
XX
ID  V$MYOD_01
XX
DT  19.10.92 (created); ewi.
DT  16.10.95 (updated); ewi.
XX
NA  MyoD
XX
DE  myoblast determination gene product
XX
BF  T00526; MyoD; Species: mouse, Mus musculus.
XX
P0      A      C      G      T
01     100     200     200     0     S
02     200     100     200     0     R
03     300     0     100     100     A
04     0     500     0     0     C
05     500     0     0     0     A
06     0     0     400     100     G
07     0     100     400     0     G
08     0     0     0     500     T
09     0     0     500     0     G
10     0     100     200     200     K
11     0     200     0     300     Y
12     100     0     300     100     G
XX
BA  5 functional elements in 3 genes
XX
//
AC  M00002
XX
ID  V$E47_01
XX
DT  19.10.92 (created); ewi.
DT  16.10.95 (updated); ewi.
XX
NA  E47
XX
DE  E47
XX
BF  T00207; E47; Species: human, Homo sapiens.
XX
P0      A      C      G      T
00     400     400     300     0     N
02     200     500     400     0     S
03     300     200     400     200     N
04     200     0     900     0     G
05     0     1100     0     0     C
06     1100     0     0     0     A
07     0     0     1100     0     G
08     100     200     800     0     G
09     0     0     0     1100     T
10     0     0     1100     0     G
11     0     0     400     700     K
12     100     400     300     300     N
13     100     600     200     200     C
14     100     400     400     200     N
15     100     400     200     300     N
XX
BA  11 selected strong binding sites for E47, E47-MyoD, E12+MyoD 
BA  and (weak) for E12
XX
CC  Group I in [903]; 5 sites selected in vitro for binding to E12N 
CC  (=N-terminally truncated E12); matrix corrected according to 
CC  the published sequences
XX
RN  [1]
RA  Sun X.-H., Baltimore D.
RT  An inhibitory domain of E12 transcription factor prevents 
RT  DNA binding in E12 homodimers but not in E12 heterodimers
RL  Cell 64:459-470 (1991).
XX
"""

def test_reader():
    input = StringIO( sample )
    motifs = list( transfac.TransfacReader( input ) )
    assert len( motifs ) == 2
    # Single value parse
    assert motifs[1].accession == "M00002"
    # Value list parse
    assert motifs[1].dates == [ '19.10.92 (created); ewi.', '16.10.95 (updated); ewi.' ]
    # Matrix parse
    assert motifs[1].matrix.sorted_alphabet == ['A','C','G','T']
    assert allclose( motifs[1].matrix.values[0], [400,400,300,0] )