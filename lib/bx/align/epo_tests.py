"tests for bx.align.epo"

import unittest, logging, pdb
import random
import numpy as np

from bx.align._epo import cummulative_intervals, bed_union
from bx.align.epo import *

class TestBed( unittest.TestCase ):
    def setUp(self):
        self.N = random.randint(1, 1000)

    def test_ci(self):
        S, D = [], []
        for i in range(self.N):
            S.append( random.randint(10, 50) )
            D.append( random.randint(10, 50) )
        D[-1] = 0
        C = cummulative_intervals(np.array(S, dtype=np.int64),
                np.array(D, dtype=np.int64))
        for i in range(self.N):
            assert C[i,1] - C[i,0] == S[i]
        for i in range(1, self.N):
            assert C[i,0] - C[i-1,1] == D[i-1], "[%d] %d != %d" % (i, C[i,0] - C[i-1,1], D[i-1])

    def test_elem_u(self):
        # back to back, so should return a single interval
        EL = []
        th = 0
        for i in range(self.N):
            size = random.randint(1, 20)
            EL.append( (th, th+size) )
            th += size
        U = bed_union( np.array(EL, dtype=np.uint64) )
        assert U[0,0] == 0 and U[0,1] == th

        # disjoint
        EL = []
        th = 0
        for i in range(self.N):
            size = random.randint(1, 20)
            EL.append( (th, th+size) )
            th += (size + 1)
        U = bed_union( np.array(EL, dtype=np.uint64) )
        for i in range(U.shape[0]):
            assert (U[i,0], U[i,1]) == EL[i]

        # random with some empty elements
        EL = []
        th = 0
        for i in range(self.N):
            size = random.randint(1, 20)
            EL.append( (th, th+size) )
            th += random.randint(1, size+size) #50% of overlapping
        U = bed_union( np.array(EL, dtype=np.uint64) )

        assert U[0,1] > U[0,0]
        for i in range(1, U.shape[0]):
            assert U[i,1] > U[i,0]
            assert U[i,0] > U[i-1,1]

cigar_pairs = [
("GGACCTGGAGAGATCAG---------------------------GACTTCAACTGTGTG-------------TCTTAGACTGGG--------AGGGTGTTA",
 "AGGCCAGGAGAGATCAGGTAAGTCTTAATTTAATAAAGAGATAGGACCTGAACTGTGTCTAACAATAGGTAATATTAGACTGGGGGAGAGAGAAGACTTTC"),

("TTT--------------------------------------------------------------------------------------------------------------------T",
 "CTTGTACCAAGGACAGTACTGGCAGCCTAATTGCTAACACTTTGTGGTGGATTGGTCCACTCAATATTTGTTCCCACCTCTTTTCAGTCCAGTTCTATAAAGGACAGAAAGTTGAAAACT"),

("A-------------------------------------------------ACACTGGACACAGCACTAACACGATTACTTA",
 "ACATTTCCCACACTCCCTTGCAGCTAGGTTTCTAGATATAATTTAGATTCCA----------------------------A"),

("TTTGGTCCTCTGGA------CGAGCAGCCAGTGCT---------------------------------------------------------------------------AAAAAAAA",
 "T---CATTCTAGCAGGTGCTGCAGCAGCAGGTAGCCCTGGAGCCAACAGTTGTGGCTATGATTCTTGATCATCAGATTTGGCTCAAGTGATGTGTTCCTCTAGCATGCACTTGAGATA"),

("G-----------------------C----------------------------------------------------------------------------------------A",
 "GGCCTGCACTGCCAGTAATTTTAACAAATTTTTAGGCACTGAATTCCCTGTATTAAATCTGTTTTCCTTAGCGTAAACAGATCTCTGTTAAATGAAACTAAACCCTGACTGATA"),

("TATT----------------------------------T",
 "TCCTTCATTTTATTTCTCCCTTAAAATTTTTTTTATTACT"),

("TAAAAA--A------A------------------------------------------------------------TTTTTTTTTTT",
 "T---AATTATTTTGCAGCAGGTCCTTGATAACATATCATCTATAAATATTTCAGCAAGAATCTCTAAAAGGCAAGAACCTCCTTCTT"),

("AAACAA---------------------------------------TT---T",
 "AAACAATACCACTGCATCACTATCAAACCCAAAAAATAACAAAAATTGGGT"),

("TCTTAAC---TGCTGAGCCATCCCTCCAGCTCCTGTTTTATTTTTATTATGAAGTAATAATA--ATAG--TAATAATAATGATG",
 "TACACTTAATTCTAAAACTTGTTATGAATCATCA----------TTGG--TTTTTTATTGTGAAGAACTAATATAATCAGA--G"),

("ATGATAATGGTATCCTAGCTCAACACCTG-GAGTTCACCCCAACAGTTAACTAA----GTTTGAGGAAGTGTTAACAAGCCTA---ACAAAGAGGACATGCCAATAGCTGACAGAGTCAC",
 "A-------CCTCTGCTAGCTCAACTCCTGAGAATCAATTATATAAGCTAGGTCAGTGGTTTTGAGAAAGTATTAGTAGACATTTCTCCAAAGAATACATAAAAATGGCC-A--CAAGTAT")
]

def toCigar(species, id, s):
    I = [(0,0)]
    L = map(len, s.split("-"))
    NZ = filter(lambda _: _, L)

    if L[0] > 0:
        I.append( (0, L[0]) )
        NZ = NZ[1:]
        L = L[1:]

    for i in range(len(NZ)):
        L.insert(0, 0)
        size = NZ[i]
        start = L.index(size)
        I.append( (I[-1][1] + start, I[-1][1]+start+size) )
        L = L[start+1:]
    if len(L):
        I.append((I[-1][1] + len(L), I[-1][1] + len(L)))
    #print I[1:]
    C = []
    for i in range(1,len(I)):
        dl = I[i][0] - I[i-1][1]
        ml = I[i][1] - I[i][0]

        dc = ""
        if dl:
            dc = (dl > 1 and str(dl) or "") + "D"

        mc = ""
        if ml:
            mc = (ml > 1 and str(ml) or "") + "M"

        C.append( dc+mc )
    MSUM = sum(map(lambda i: i[1]-i[0], I))
    start = random.randint(50, 10000)
    return "%s\t%d\t1\t%d\t%d\t%d\t%s" % (species, id, start, start+MSUM-1,
            random.choice((-1,1)), "".join(C))



class TestEpo( unittest.TestCase ):
    def setUp(self):
        self.epo_records = []
        for i, (t,q) in enumerate(cigar_pairs):
            gab_pair = (toCigar("homo_sapiens", i, t), toCigar("mus_musculus", i, q))
            A = EPOitem._strfactory(gab_pair[0])
            B = EPOitem._strfactory(gab_pair[1])
            if A and B:
                self.epo_records.append( (A,B) )

    def test_out(self):
        def ch(c, ci):
            th = 0
            for l,t in ci:
                if t == 'M':
                    assert c[th:th+l].find('-') == -1
                else:
                    assert c[th:th+l] == '-' * l
                th += l

        for (a,b) in self.epo_records:
            ca, cb = cigar_pairs[int(a.gabid)]
            #if a.strand == '-': ca = ca[::-1]
            #if b.strand == '-': cb = cb[::-1]
            ch(ca, a.cigar_iter(False))
            ch(cb, b.cigar_iter(False))

    def test_make_chain(self):
        def cch(cigar, s, e):
            return cigar[s:e].find('-') == -1

        for p in self.epo_records:
            chain  = Chain._make_from_epo(p[0], p[1], {"chr1" : 500}, {"chr1" : 800})
            if not chain:
                continue
            ch, S, T, Q = chain
            i = int( ch.id )
            c1, c2 = cigar_pairs[i]
            if p[0].strand == '-':
                c1 = c1[::-1]
                c2 = c2[::-1]
            th = 0
            for s, t, q in zip(S, T, Q):
                if not (cch(c1, th, th+s) and cch(c2, th, th+s)):
                    pdb.set_trace()
                assert cch(c1, th, th+s) and cch(c2, th, th+s), "%s and %s" % (c1[th:th+s], c2[th:th+s])
                if t > q:
                    cch(c1, th+s, th+s+t) and c1[th+s:th+s+t] == '-'*t
                else:
                    cch(c2, th+s, th+s+q) and c1[th+s:th+s+q] == '-'*q
                th = th + s + max(t,q)

    def test_rem_dash(self):
        # ****--****-------****  4M2D4M7D4M
        # *******-------*******  7M7D7M
        # has 4 dash columns and should become
        # ****--****---****      4M2D4M3D4M
        # *******---*******      7M3D7M

        for i in range(100):
            dash_cols = random.randint(0, 10)
            tStart = random.randint(0, 1000)
            qStart = random.randint(0, 1000)
            epo_pair = (EPOitem._strfactory("homo_sapiens\t0\t1\t%d\t%d\t1\t%s" % (tStart, tStart+12-1, "4M2D4M%dD4M" % (dash_cols+3))),
                    EPOitem._strfactory("mus_musculus\t0\t1\t%d\t%d\t1\t%s" % (qStart, qStart+14-1, "7M%dD7M" % (dash_cols+3))))
            chain = Chain._make_from_epo(epo_pair[0], epo_pair[1], {"chr1":500}, {"chr1":800})
            ti = epo_pair[0].intervals(False)
            qi = epo_pair[1].intervals(False)
            assert ti[2][0] - ti[1][1] - dash_cols == chain[2][1]
            assert qi[1][0] - qi[0][1] - dash_cols == chain[2][1]

        # ----*****
        # *-------*
        # has 3 dash cols and should become
        # *
        # *
        # with the qStart += 1 and tStart += 4

        for i in range(100):
            dash_cols = random.randint(0, 10)
            tm = random.randint(6, 10)
            qm = random.randint(1, 5)

            tStart = random.randint(0, 1000)
            qStart = random.randint(0, 1000)

            epo_pair = (EPOitem._strfactory("homo_sapiens\t0\t1\t%d\t%d\t1\t%s" % (tStart, tStart+tm-1,
                    "%dD%dM" % (dash_cols+1, tm))),
                EPOitem._strfactory("mus_musculus\t0\t1\t%d\t%d\t1\t%s" % (qStart, qStart+qm+1-1,
                    "M%dD%dM" % (dash_cols+tm-qm, qm))))
            chain = Chain._make_from_epo(epo_pair[0], epo_pair[1], {"chr1":500}, {"chr1":800})
            if chain[1][-1] != qm:
                pdb.set_trace()
            assert chain[1][-1] == qm
            # correct also for coordinate interpretation differences between UCSC and EPO
            assert (qStart + 1) - 1 == chain[0].qStart, "%d != %d" % (qStart + 1, chain[0].qStart)


if __name__ == '__main__':
    unittest.main()
