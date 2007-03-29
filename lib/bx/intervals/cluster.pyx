"""
ClusterTree data structure that supports efficient queries for "cluster" of
intervals (groups where no interval is more than some distance from at least
one other).

TODO: Documentation of the algorithmic details.
"""

cdef extern from "common.h":
    pass

cdef extern from "cluster.h":
    cdef struct listitem:
        int value
        listitem * next

    cdef struct linelist:
        listitem * head
        listitem * tail
        
    cdef struct ClusterNode:
        int start
        int end
        int regions
        ClusterNode * left
        ClusterNode * right
        linelist * linenums

    cdef struct treeitr:
        ClusterNode * next
        ClusterNode * value

    # Insert into tree
    ClusterNode * clusterNodeInsert(ClusterNode ** cn, int start, int end, int linenum, int mincols)
    
    # Find node holding position
    ClusterNode * clusterNodeFind(ClusterNode * cn, int position)
    # Free entire tree
    void freetree(ClusterNode ** cn)

    void get_itr_in(ClusterNode * root, treeitr **itr)

    ClusterNode * next(treeitr **itr)

    int has_next(treeitr **itr)
    
cdef class ClusterTree:
    cdef ClusterNode * root
    cdef int mincols
    cdef int minregions
    
    def __new__( self, mincols, minregions ):
        self.root = NULL
        self.mincols = mincols
        self.minregions = minregions

    def __dealloc__( self ):
        freetree( &self.root )
        
    def insert( self, int start, int end, int linenum):
        clusterNodeInsert( &self.root, start, end, linenum, self.mincols )
        return self

    def getregions( self ):
        cdef treeitr * myitr
        cdef ClusterNode * node
        cdef linelist* nums
        cdef listitem* num

        myitr = NULL
        get_itr_in(self.root, &myitr)
        regions = []
        while (has_next(&myitr)):
            node = next(&myitr)
            if node.regions >= self.minregions:
                lines = []
                if node.linenums and node.linenums.head:
                    nums = node.linenums
                    num = nums.head
                    while num:
                        lines.append( num.value )
                        num = num.next
                regions.append( (node.start, node.end, lines) )

        return regions
    
    def getlines( self ):
        cdef treeitr * myitr
        cdef ClusterNode * node
        cdef linelist* nums
        cdef listitem* num
        
        myitr = NULL
        get_itr_in(self.root, &myitr)
        lines = []
        while (has_next(&myitr)):
            node = next(&myitr)
            if node.regions >= self.minregions:
                if node.linenums and node.linenums.head:
                    nums = node.linenums
                    num = nums.head
                    while num:
                        lines.append( num.value )
                        num = num.next
        return lines
