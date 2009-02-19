"""
Data structure for performing intersect queries on a set of intervals which
preserves all information about the intervals (unlike bitset projection methods).

:Authors: James Taylor (james@jamestaylor.org),
          Ian Schenk (ian.schenck@gmail.com),
          Brent Pederson (bpederse@gmail.com)
"""

# Historical note:
#    This module original contained an implementation based on sorted endpoints
#    and a binary search, using an idea from Scott Schwartz and Piotr Berman.
#    Later an interval tree implementation was implemented by Ian for Galaxy's
#    join tool (see `bx.intervals.operations.quicksect.py`). This was then
#    converted to Cython by Brent, who also added support for
#    upstream/downstream/neighbor queries. This was modified by James to
#    handle half-open intervals strictly, to maintain sort order, and to
#    implement the same interface as the origianl Intersecter.

import operator

cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

ctypedef char * char_star

cdef class Interval:
    """
    Basic feature, with required integer start and end properties.
    Also accpets optional strand as +1 or -1 (used for up/downstream queries),
    a name, and any arbitrary data is sent in on the info keyword argument

    >>> from bx.intervals.intersection import Interval

    >>> f1 = Interval(23, 36)
    >>> f2 = Interval(34, 48, strand=-1, name="fred", info={'chr':12, 'anno':'transposon'})
    >>> f2
    Interval(34, 48, strand=-1, name="fred", {'anno': 'transposon', 'chr': 12})

    """
    cdef public int start, end, strand
    cdef public object info
    cdef public object name

    def __init__(self, int start, int end, int strand=0, int chr=0, object name="", object info=None):
        assert start <= end, "start must be less than end"
        self.start  = start
        self.end   = end      
        self.strand = strand
        self.name   = name
        self.info   = info

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if self.strand != 0:
            fstr += ", strand=%d" % self.strand
        if strlen(self.name) != 0:
            fstr += ', name="' + str(self.name) + '"'
        if not self.info is None:
            fstr += ", " + str(self.info)
        fstr += ")"
        return fstr

    def __cmp__(self, other):
        return cmp( self.start, other.start ) or cmp( self.end, other.end )

cdef inline int imax2(int a, int b):
    if b > a: return b
    return a

cdef inline int imax3(int a, int b, int c):
    if b > a: 
        if c > b:
            return c
        return b
    if a > c:
        return a
    return c

cdef inline int imin3(int a, int b, int c):
    if b < a: 
        if c < b:
            return c
        return b
    if a < c:
        return a
    return c

cdef inline int imin2(int a, int b):
    if b < a: return b
    return a

cdef float nlog = -1.0 / log(0.5)

cdef class IntervalNode:
    """
    Data structure for performing intersect and neighbor queries on a 
    set of intervals. Algorithm uses a segment/interval tree to perform
    efficient queries. 

    Usage
    =====
    >>> from bx.intervals.intersection import IntervalNode, Interval
    >>> tree = IntervalNode(Interval(0, 10, -1))

    Add intervals, the only requirement is that the interval have integer
    start and end attributes. Optional arguments are strand, name, and info.

    >>> Interval(1, 22, strand=-1, name="fred", info={'chr':12, 'anno': 'anything'})
    Interval(1, 22, strand=-1, name="fred", {'anno': 'anything', 'chr': 12})


    >>> tree = tree.insert(Interval(3, 7, 1))
    >>> tree = tree.insert(Interval(3, 40, -1))
    >>> tree = tree.insert(Interval(13, 50, 1))

    Queries
    -------

    find
    ++++

    >>> tree.find(2, 5)
    [Interval(0, 10, strand=-1), Interval(3, 7, strand=1), Interval(3, 40, strand=-1)]
    >>> tree.find(7, 5)
    [Interval(0, 10, strand=-1), Interval(3, 40, strand=-1)]
    >>> tree.find(11, 100)
    [Interval(3, 40, strand=-1), Interval(13, 50, strand=1)]
    >>> tree.find(100, 200)
    []

    left/right
    ++++++++++
    the left method finds features that are strictly to the left of
    the query feature. overlapping features are not considered:

    >>> tree.left(Interval(0, 1))
    []
    >>> tree.left(Interval(11, 12))
    [Interval(0, 10, strand=-1)]


    up/downstream
    +++++++++++++
    up/downstream method behave exactly like left/right, except that
    the direction is determined by the strand of the query feature. 
    If the strand is 1, then upstream is left, downstream is right.

    If the strand is -1, then upstream is right, downstream is left.
    >>> tree.upstream(Interval(11, 12, strand=1))
    [Interval(0, 10, strand=-1)]
    >>> tree.upstream(Interval(11, 12, strand=-1))
    [Interval(13, 50, strand=1)]

    all of these method take an argument 'n' for the number of results desired.
    >>> tree.upstream(Interval(1, 2, strand=-1), n=3)
    [Interval(3, 7, strand=1), Interval(3, 40, strand=-1), Interval(13, 50, strand=1)]
    
    nearest neighbors
    +++++++++++++++++
    #>>> tree.nearest_neighbors(Interval(1, 2))
    #[Interval(0, 10, strand=-1)]

    #>>> tree.nearest_neighbors(Interval(1, 2), n=2)
    #[Interval(0, 10, strand=-1), Interval(3, 7, strand=1)]

    """
    cdef float priority
    cdef Interval interval 
    cdef public int start, end
    cdef int minend, maxend, minstart
    cdef IntervalNode cleft, cright, croot

    property left_node:
        def __get__(self):
            return self.cleft if self.cleft is not EmptyNode else None
    property right_node:
        def __get__(self):
            return self.cright if self.cright is not EmptyNode else None
    property root_node:
        def __get__(self):
            return self.croot if self.croot is not EmptyNode else None
    
    def __repr__(self):
        return "IntervalNode(%i, %i)" % (self.start, self.end)

    def __cinit__(IntervalNode self, Interval interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority   = ceil(nlog * log(-1.0/(1.0 * rand()/RAND_MAX - 1)))
        self.start      = interval.start
        self.end       = interval.end
        self.interval   = interval
        self.maxend    = interval.end
        self.minstart   = interval.start
        self.minend    = interval.end
        self.cleft       = EmptyNode
        self.cright      = EmptyNode
        self.croot       = EmptyNode

    def insert(self, interval):
        return self._insert(interval)

    cdef IntervalNode _insert(IntervalNode self, Interval interval):
        cdef IntervalNode croot = self
        # If starts are the same, decide which to add interval to based on
        # end, thus maintaining sortedness relative to start/end
        cdef int decision_endpoint = interval.start
        if interval.start == self.start:
            decision_endpoint = interval.end
        
        if decision_endpoint > self.start:
            # insert to cright tree
            if self.cright is not EmptyNode:
                self.cright = self.cright._insert(interval )
            else:
                self.cright = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.cright.priority:
                croot = self.rotate_left()
        else:
            # insert to cleft tree
            if self.cleft is not EmptyNode:
                self.cleft = self.cleft._insert(interval)
            else:
                self.cleft = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.cleft.priority:
                croot = self.rotate_right()
    
        croot.set_ends()
        self.cleft.croot  = croot
        self.cright.croot = croot
        return croot

    cdef IntervalNode rotate_right(IntervalNode self):
        cdef IntervalNode croot = self.cleft
        self.cleft  = self.cleft.cright
        croot.cright = self
        self.set_ends()
        return croot

    cdef IntervalNode rotate_left(IntervalNode self):
        cdef IntervalNode croot = self.cright
        self.cright = self.cright.cleft
        croot.cleft  = self
        self.set_ends()
        return croot

    cdef inline void set_ends(IntervalNode self):
        if self.cright is not EmptyNode and self.cleft is not EmptyNode: 
            self.maxend = imax3(self.end, self.cright.maxend, self.cleft.maxend)
            self.minend = imin3(self.end, self.cright.minend, self.cleft.minend)
            self.minstart = imin3(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxend = imax2(self.end, self.cright.maxend)
            self.minend = imin2(self.end, self.cright.minend)
            self.minstart = imin2(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxend = imax2(self.end, self.cleft.maxend)
            self.minend = imin2(self.end, self.cleft.minend)
            self.minstart = imin2(self.start, self.cleft.minstart)
        

    def intersect(self, int start, int end, sort=True ):
        """
        given a start and a end, return a list of features
        falling within that range
        """
        cdef list results = []
        self._intersect(start, end, results)
        return results

    find = intersect
        
    cdef void _intersect(IntervalNode self, int start, int end, list results):
        # Left subtree
        if self.cleft is not EmptyNode and self.cleft.maxend > start:
            self.cleft._intersect(start, end, results)
        # This interval
        if (self.end > start) and (self.start < end): results.append(self.interval)
        # Right subtree
        if self.cright is not EmptyNode and self.start < end:
            self.cright._intersect(start, end, results)
    

    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend + max_dist < position: return
        if self.minstart > position: return

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cright is not EmptyNode:
                self.cright._seek_left(position, results, n, max_dist)

        if -1 < position - self.end < max_dist:
            results.append(self.interval)

        # TODO: can these conditionals be more stringent?
        if self.cleft is not EmptyNode:
                self.cleft._seek_left(position, results, n, max_dist)


    
    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend < position: return
        if self.minstart - max_dist > position: return

        #print "SEEK_RIGHT:",self, self.cleft, self.maxend, self.minstart, position

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cleft is not EmptyNode: 
                self.cleft._seek_right(position, results, n, max_dist)

        if -1 < self.start - position < max_dist:
            results.append(self.interval)

        if self.cright is not EmptyNode:
                self.cright._seek_right(position, results, n, max_dist)

    def neighbors(self, Interval f, int n=1, int max_dist=2500):
        cdef list neighbors = []

        cdef IntervalNode right = self.cright
        while right.cleft is not EmptyNode:
            right = right.cleft

        cdef IntervalNode left = self.cleft
        while left.cright is not EmptyNode:
            left = left.cright
        return [left, right]

    cpdef left(self, Interval f, int n=1, int max_dist=2500):
        """find n features with a start > than f.end
        f: a Interval object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use start - 1 becuase .left() assumes strictly left-of
        self._seek_left(f.start - 1, results, n, max_dist)
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('end'), reverse=True)
        return r[:n]

    cpdef right(self, Interval f, int n=1, int max_dist=2500):
        """find n features with a end < than f.start
        f: a Interval object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use end + 1 becuase .right() assumes strictly right-of
        self._seek_right(f.end + 1, results, n, max_dist)
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        return r[:n]

    def upstream(self, Interval f, int n=1, int max_dist=2500):
        """find n upstream features where upstream is determined by
        the strand of the query Interval f
        Overlapping features are not considered.

        f: a Interval object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        if f.strand == -1:
            return self.right(f, n, max_dist)
        return self.left(f, n, max_dist)

    def downstream(self, Interval f, int n=1, int max_dist=2500):
        """find n downstream features where downstream is determined by
        the strand of the query Interval f
        Overlapping features are not considered.

        f: a Interval object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        if f.strand == -1:
            return self.left(f, n, max_dist)
        return self.right(f, n, max_dist)

    def traverse(self, func):
        self._traverse(func)

    cdef void _traverse(IntervalNode self, object func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self)
        if self.cright is not EmptyNode: self.cright._traverse(func)

cdef IntervalNode EmptyNode = IntervalNode(Interval(0, 0))

## ---- Wrappers that retain the old interface -------------------------------

cdef class WrapperInterval( Interval ):
    pass

cdef class Intersecter:
    cdef IntervalNode intervals
    """
    Data structure for performing window intersect queries on a set of 
    intervals. Now a wrapper to IntervalNode

    Usage
    =====

    >>> from bx.intervals.intersection import Intersecter, Interval
    >>> intersecter = Intersecter()

    Add intervals, the only requirement is that the interval have integer
    start and end attributes:

    >>> intersecter.add_interval( Interval( 0,  10 ) )
    >>> intersecter.add_interval( Interval( 3,  7 ) )
    >>> intersecter.add_interval( Interval( 3,  40 ) )
    >>> intersecter.add_interval( Interval( 10, 50 ) )

    Perform queries:

    >>> intersecter.find( 2, 5 )
    [Interval(0, 10), Interval(3, 7), Interval(3, 40)]
    >>> intersecter.find( 10, 100 )
    [Interval(3, 40), Interval(10, 50)]
    >>> intersecter.find( 100, 200 )
    []
    
    """

    def __cinit__( self ):
        """Initialize"""
        self.intervals = None

    cpdef add_interval( self, interval ):
        """Add an interval to the stored set"""
        ## assert interval.start < interval.end, "Intervals must have length >= 1"
        # If not given a feature, we must make one:
        if not isinstance( interval, Interval ):
            interval = WrapperInterval( interval.start, interval.end, info=interval )
        if self.intervals is None:
            self.intervals = IntervalNode( interval )
        else:
            self.intervals = self.intervals.insert( interval )

    cpdef add( self, start, end, value=None ):    
        self.add_interval( Interval( start, end, value ) )
        
    cpdef find( self, start, end ):
        """Return a list of all stored intervals intersecting [start,end)"""
        rval = self.intervals.find( start, end )
        for i in range( len( rval ) ):
            if isinstance( rval[i], WrapperInterval ):
                rval[i] = rval[i].info
        rval.sort()
        return rval