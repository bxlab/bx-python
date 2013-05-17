
"""Classes and utilities for mutliple alignments from the EPO pipeline"""

from __future__ import with_statement

import logging, re, cPickle, os
from collections import namedtuple
from itertools import imap

from _epo import rem_dash, fastLoadChain, bed_union, cummulative_intervals


log = logging.getLogger(__name__)

class Chain( namedtuple('Chain', 'score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id') ):
    """A Chain header as in http://genome.ucsc.edu/goldenPath/help/chain.html

    chain coordinates are with respect to the strand, so for example tStart on the + strand is the
    distance from the leftmost position; tStart on the - strand is the distance from the rightmost position."""

    __slots__ = ()

    def __str__(self):
        return "chain {score} {tName} {tSize} {tStrand} {tStart} {tEnd} {qName} {qSize} {qStrand} {qStart} {qEnd} {id}".format(**self._asdict())

    @classmethod
    def _strfactory(cls, line):
        """factory class method for Chain

        :param line: header of a chain (in .chain format)
        """

        assert type(line) == str, "this is a factory from string"

        line = line.rstrip().split()[1:] # the first component is the keyword "chain"
        tup = map(lambda t: t[0](t[1]),
                zip([int, str, int, str, int, int, str, int, str, int, int, str], line))
        return tuple.__new__(cls, tup)

    @classmethod
    def _make_from_epo(cls, trg_comp, qr_comp, trg_chrom_sizes, qr_chrom_sizes):
        """crate a chain of collinear rings from the given  components.

        The target of the chain will always be on the forward strand.
        This is done to avoid confusion when mapping psl files. So,
        if trg_comp.strand=-, qr_comp.strand=- (resp. +) the
        chain header will have tStrand=+, qStrand=+ (resp. -). No strand
        changes on the other cases.

        :param trg_comp: target (i.e, the first) component
        :type trg_comp: L{EPOitem}
        :param qr_comp: query (i.e, the second) component
        :type qr_comp: L{EPOitem}
        :param trg_chrom_sizes: chromosome sizes of the target
        :type trg_chrom_sizes: dictionary of the type (chrom) --> size
        :param qr_chrom_sizes: chromosome sizes of the query
        :type qr_chrom_sizes: dictionary of the type (chrom) --> size
        :return: A L{Chain} instance"""

        # size, target, query arrays
        S, T, Q = [], [], []

        #the target strand of the chain must be on the forward strand
        trg_intervals = trg_comp.intervals(reverse = trg_comp.strand == '-')
        qr_intervals = qr_comp.intervals(reverse = trg_comp.strand == '-')
        if len(trg_intervals) == 0 or len(qr_intervals) == 0:
            log.warning("deletion/insertion only intervals")
            return None
        A, B = rem_dash(trg_intervals, qr_intervals)
        # correct for when cigar starts/ends with dashes (in number of bases)
        tr_start_correction = max(B[0][0] - A[0][0], 0)
        tr_end_correction = max(A[-1][1] - B[-1][1], 0)
        qr_start_correction = max(A[0][0] - B[0][0], 0)
        qr_end_correction = max(B[-1][1] - A[-1][1], 0)

        a, b = A.pop(0), B.pop(0)

        # intervals are 0-base, halfo-open => lengths = coordinate difference
        while A or B:
            if a[1] < b[1]:
                T.append(0); Q.append( A[0][0] - a[1] ); S.append( min(a[1], b[1]) - max(a[0], b[0]) )
                a = A.pop(0)
            elif b[1] < a[1]:
                Q.append(0); T.append( B[0][0] - b[1] ); S.append( min(a[1], b[1]) - max(a[0], b[0]) )
                b = B.pop(0)
            elif A and B:
                assert 1 > 2, "there are dash columns"
            else:
                break
        S.append( min(a[1], b[1]) - max(a[0], b[0]) )
        assert len(T) == len(Q) == len(S) - 1, "(S, T, Q) = (%d, %d, %d)" % tuple(map(len, (S, T, Q)))

        tSize = trg_chrom_sizes[trg_comp.chrom]
        qSize = qr_chrom_sizes[qr_comp.chrom]
        ## UCSC coordinates are 0-based, half-open and e! coordinates are 1-base, closed
        ## chain_start = epo_start - 1 and chain_end = epo_end
        if qr_comp.strand == '+':
            chain = Chain(0,
                    trg_comp.chrom, tSize, "+",
                    (trg_comp.start - 1) + tr_start_correction, trg_comp.end  - tr_end_correction,
                    qr_comp.chrom, qSize, (qr_comp.strand == trg_comp.strand and '+' or '-'),
                    (qr_comp.start - 1) + qr_start_correction, qr_comp.end  - qr_end_correction,
                    qr_comp.gabid)
        else:
            chain = Chain(0,
                    trg_comp.chrom, tSize, "+",
                    (trg_comp.start - 1) + tr_start_correction, trg_comp.end - tr_end_correction,
                    qr_comp.chrom, qSize, (qr_comp.strand == trg_comp.strand and '+' or '-'),
                    (qr_comp.start - 1) + qr_end_correction, qr_comp.end - qr_start_correction,
                    qr_comp.gabid)

        # strand correction. in UCSC coordinates this is: size - coord
        if chain.qStrand == '-':
            chain = chain._replace(qEnd = chain.qSize - chain.qStart,
                    qStart = chain.qSize - chain.qEnd)

        assert chain.tEnd - chain.tStart  == sum(S) + sum(T), "[%s] %d != %d" % (str(chain),
                chain.tEnd - chain.tStart, sum(S) + sum(T))
        assert chain.qEnd - chain.qStart == sum(S) + sum(Q), "[%s] %d != %d" % (str(chain),
                chain.qEnd - chain.qStart, sum(S) + sum(Q))
        return chain, S, T, Q

    def slice(self, who):
        "return the slice entry (in a bed6 format), AS IS in the chain header"

        assert who in ('t', 'q'), "who should be 't' or 'q'"

        if who == 't':
            return (self.tName, self.tStart, self.tEnd, self.id, self.score, self.tStrand)
        else:
            return (self.qName, self.qStart, self.qEnd, self.id, self.score, self.qStrand)

    def bedInterval(self, who):
        "return a BED6 entry, thus DOES coordinate conversion for minus strands"

        if who == 't':
            st, en = self.tStart, self.tEnd
            if self.tStrand == '-':
                st, en = self.tSize-en, self.tSize-st
            return (self.tName, st, en, self.id, self.score, self.tStrand)
        else:
            st, en = self.qStart, self.qEnd
            if self.qStrand == '-':
                st, en = self.qSize-en, self.qSize-st
                assert en-st == self.qEnd - self.qStart
            return (self.qName, st, en, self.id, self.score, self.qStrand)

    @classmethod
    def _parse_file(cls, path, pickle=False):
        """parse a .chain file into a list of the type [(L{Chain}, arr, arr, arr) ...]

        :param fname: name of the file"""

        fname = path
        if fname.endswith(".gz"):
            fname = path[:-3]

        if fname.endswith('.pkl'):
            #you asked for the pickled file. I'll give it to you
            log.debug("loading pickled file %s ..." % fname)
            return cPickle.load( open(fname) )
        elif os.path.isfile("%s.pkl" % fname):
            #there is a cached version I can give to you
            log.info("loading pickled file %s.pkl ..." % fname)
            if os.stat(path).st_mtime > os.stat("%s.pkl" % fname).st_mtime:
                log.critical("*** pickled file %s.pkl is not up to date ***" % (path))
            return cPickle.load( open("%s.pkl" % fname) )

        data = fastLoadChain(path, cls._strfactory)
        if pickle and not os.path.isfile('%s.pkl' % fname):
            log.info("pckling to %s.pkl" % (fname))
            with open('%s.pkl' % fname, 'wb') as fd:
                cPickle.dump(data, fd)
        return data

class EPOitem(namedtuple('Epo_item', 'species gabid chrom start end strand cigar')):
    "this format is how alignments are delivered from e!"

    __slots__ = ()

    cigar_pattern = re.compile("(\d*)([MD])")

    def __repr__(self): return str(self)

    def __str__(self):
        c = self.cigar[:5] + "..." + self.cigar[-5:]
        return "(%s %s %s %d %d %s %s)" % tuple(self[:6] + (c,))

    @classmethod
    def _strfactory(cls, line):
        """factory method for an EPOitem

        :param line: a line of input"""

        cmp = line.rstrip().split()
        chrom = cmp[2]
        if not chrom.startswith("chr"):
            chrom = "chr%s" % chrom
        instance =  tuple.__new__(cls,
                (cmp[0], cmp[1],
                    chrom, int(cmp[3]), int(cmp[4]),
                    {'1' : '+', '-1' : '-'}[cmp[5]], cmp[6]))
        span = instance.end - instance.start + 1
        m_num = sum( map(lambda t: (t[1] == "M" and [t[0]] or [0])[0] , instance.cigar_iter(False)) )
        if span != m_num:
            log.warning("[{gabid}] {species}.{chrom}:{start}-{end}.".format(**instance._asdict()) + "(span) %d != %d (matches)" % (span, m_num))
            return None
        return instance

    @classmethod
    def _parse_epo(cls, fname):
        """Load an entire file in the EPO format into a dictionary of the type {gab_id => [Epoitem, ...]}

        :param fname: file name"""

        data = {}
        with open(fname) as fd:
            for el in imap(cls._strfactory, fd):
                if el:
                    data.setdefault(el.gabid, []).append( el )
        log.info("parsed %d elements from %s" % (len(data), fname))
        return data

    def cigar_iter(self, reverse):
        """self.cigar => [(length, type) ... ] iterate the cigar

        :param reverse: whether to iterate in the reverse direction (right-to-left)
        :type reverse: boolean

        :return a list of pairs of the type [(length, M/D) ..]
        """

        l = 0
        P = self.cigar_pattern

        data = []
        cigar = self.cigar
        parsed_cigar = re.findall(P, cigar)
        if reverse:
            parsed_cigar = parsed_cigar[::-1]
        for _l, t in parsed_cigar:
            # 1M is encoded as M
            l = (_l and int(_l) or 1) # int(_l) cannot be 0
            data.append( (l, t) )
        return data

    def intervals(self, reverse, thr=0):
        """return a list of (0-based half-open) intervals representing the match regions of the cigar

        for example 4MD4M2DM with reverse=False will produce [(0,4), (5,9), (11,12)]
        4MD4M2DM with reverse=True will produce [(0,1), (3,7), (8,12)] (= 12 - previous interval)

        :param reverse: whether to iterate in the reverse direction (right-to-left) (this is passed as is to self.cigar_iter)
        :type reverse: boolean
        :param thr: shift all intervals by this much
        :type thr: integer

        :return: list of pairs"""

        d = [(thr,thr)]
        dl = 0
        for tup in self.cigar_iter(reverse):
            if tup[1] == "D":
                dl = tup[0]
            else:
                s = d[-1][1] + dl
                d.append( (s, s+tup[0]) )

        assert d[0] == (thr, thr)
        # assert that nr. of Ms in the interval == sum of produced intervals
        assert sum( map(lambda t: t[0], filter(lambda t: t[1] == "M", self.cigar_iter(False))) ) == sum( map(lambda t: t[1]-t[0], d) )

        d_sum = sum( map(lambda t: t[1]-t[0], d) )
        assert self.end - self.start + 1 == d_sum, "[ (%d, %d) = %d ] != %d" % (self.start, self.end,
                self.end-self.start+1, d_sum)
        return d[1:] #clip the (thr, thr) entry


