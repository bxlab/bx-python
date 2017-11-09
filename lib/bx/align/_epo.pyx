

import logging, gzip
from collections import namedtuple
import numpy
cimport numpy

log = logging.getLogger(__name__)

cimport cython

DTYPE = numpy.uint64


cdef inline int max2( int a, int b ):
    if b > a:
        return b
    return a

cdef inline int min2( int a, int b ):
    if b < a:
        return b
    return a


def rem_dash(p, q):
    """remove dash columns and shift match intervals to the left. both iterables
    are read on the same direction left-to-right.
    """

    def myp(l):
        if l: return l.pop(0)

    def adv(queue, i, d):
        # shifted interval
        shi = i[0]-d, i[1]-d
        assert shi[0] >= 0
        if queue and queue[-1][1] == shi[0]:
            # join to the preceeding one
            queue[-1] = (queue[-1][0], shi[1])
        else:
            queue.append( shi )
        return queue

    p_card = sum( map(lambda i: p[i][1] - p[i][0], range(len(p))) )
    q_card = sum( map(lambda i: q[i][1] - q[i][0], range(len(q))) )

    P, Q = [], []
    dash = 0 # dash (on both cigars) count so far
    a, b = p.pop(0), q.pop(0)
    #while p or q:
    while a and b:
        assert dash <= min(a[0], b[0])
        i = max(a[0], b[0]) - min(a[1], b[1])
        if i >= 0: # no intersection
            if a[1] <= b[0]:
                if p:
                    i = min(i, p[0][0] - a[1])
                P = adv(P, a, dash)
                a = myp(p)
            else:
                if q:
                    i = min(i, q[0][0] - b[1])
                Q = adv(Q, b, dash)
                b = myp(q)
            dash += i
        else: # intersection
            if a[1] >= b[1]:
                Q = adv(Q, b, dash); b = myp(q)
            elif a[1] < b[1]:
                P = adv(P, a, dash); a = myp(p)
        #if not a or not b: # no more matchings
        #    break
    assert (not p) or (not q), "one or both should be empty: p=%s, q=%s" % (str(p), str(q))

    if a: P = adv(P, a, dash)
    if b: Q = adv(Q, b, dash)

    # remaining intervals (in q or p)
    r, R = p, P
    if q: r, R = q, Q
    # just extend the last inteval by the remaining bases
    R[-1] = (R[-1][0], R[-1][1] + sum( map(lambda i: i[1]-i[0], r) ))

    P_card = sum( map(lambda i: P[i][1] - P[i][0], range(len(P))) )
    Q_card = sum( map(lambda i: Q[i][1] - Q[i][0], range(len(Q))) )

    assert p_card == P_card, "%d != %d" % (p_card, P_card)
    assert q_card == Q_card, "%d != %d" % (q_card, Q_card)

    return P, Q

def fastLoadChain(fname, hf):
    data = []
    open_f = (fname.endswith(".gz") and gzip.open or open)
    with open_f(fname, "rt") as fd:
        while True:
            line = fd.readline()
            if line == "":
                break
            hd = hf(line)
            N = []
            line = fd.readline().split()
            while len(line) == 3:
                N.append( (int(line[0]), int(line[1]), int(line[2])) )
                line = fd.readline().split()
            if len(line) != 1:
                raise ValueError("last matching block expected (found %s)" % str(line))
            N.append( (int(line[0]), 0, 0) )
            s, t, q = zip( *N )
            data.append( (hd,
                numpy.array(s, dtype=numpy.int),
                numpy.array(t, dtype=numpy.int),
                numpy.array(q, dtype=numpy.int)) )
            assert hd.tEnd - hd.tStart == sum(s) + sum(t)
            assert hd.qEnd - hd.qStart == sum(s) + sum(q)
            fd.readline() # a blank line
        log.info("parsed %d elements from %s" % (len(data), fname))
    return data



@cython.wraparound(False)
@cython.boundscheck(False)
cpdef numpy.ndarray[numpy.uint64_t, ndim=2] bed_union( numpy.ndarray[numpy.uint64_t, ndim=2] elements ):
    """compute the union of these elements. simply walk the sorted elements and join the intersecting ones
    works on half-open intervals, i.e., [a, b), [b, c) ---> [a, c)

    @param elements: 2-dim numpy array of unsigned64 ints
    @return: 2-dim numpy array of unsigned64 ints"""

    assert numpy.shape(elements)[0] > 0

    cdef Py_ssize_t cst, cen, i, j
    cdef numpy.ndarray[numpy.uint64_t, ndim=2] tmp_elems, final_elems

    elements.sort(axis=0)
    assert elements[0][0] <= elements[numpy.shape(elements)[0]-1][0]
    tmp_elems = numpy.zeros((numpy.shape(elements)[0], 2), dtype=DTYPE)
    cst = elements[0, 0]
    cen = elements[0, 1]
    j = 0

    for i in range(1, numpy.shape(elements)[0]):
        if elements[i, 0] <= cen: # overlaps with the last one
            cen = max2(cen, elements[i, 1])
        else:
            tmp_elems[j, 0] = cst
            tmp_elems[j, 1] = cen
            j += 1
            cst = elements[i, 0]
            cen = elements[i, 1]
    tmp_elems[j, 0] = cst
    tmp_elems[j, 1] = cen
    j += 1
    final_elems = numpy.empty((j, 2), dtype=DTYPE)
    for i in range(j):
        final_elems[i, 0] = tmp_elems[i, 0]
        final_elems[i, 1] = tmp_elems[i, 1]
    assert final_elems[0, 0] == elements[0, 0], "fe=%d, e=%d" % (final_elems[0,0], elements[0,0])
    return final_elems

#@cython.wraparound(False)
#@cython.boundscheck(False)
cpdef numpy.ndarray[numpy.int64_t, ndim=2] cummulative_intervals(numpy.ndarray[numpy.int64_t, ndim=1] S,
        numpy.ndarray[numpy.int64_t, ndim=1] D ):
    """compute cummulative intervals for this side of an aligmnent. S and D are one side of
    the alignment as described in the chain file format"""

    cdef int N = S.shape[0]
    cdef int i = 0, j = 0
    assert N  == D.shape[0]
    cdef numpy.ndarray[numpy.int64_t, ndim=2] cumm_i = numpy.empty((N, 2), dtype=numpy.int64)

    cumm_i[0,0] = 0
    cumm_i[0,1] = S[0]
    for i in range(N-1):
        j = i + 1
        cumm_i[j,0] = cumm_i[i, 1] + D[i]
        cumm_i[j,1] = cumm_i[j,0] + S[j]
    return cumm_i



