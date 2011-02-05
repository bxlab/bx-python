"""
Kanwei Li, 2009
Inspired by previous ClusterTree

Provides a ClusterTree data structure that supports efficient finding of 
clusters of intervals that are within a certain distance apart.

This clustering algorithm uses a binary tree structure. Nodes correspond to 
non-overlapping intervals, where overlapping means that the distance between
two intervals is less or equal to the max separation.

The tree self-balances using rotations based on the binomial sequence. Merges
among nodes are performed whenever a node is changed/added that will cause other
nodes to form a new cluster.

C source code is in src/cluster.c
"""

cdef extern from "cluster.h":
    
    cdef struct struct_interval:
        int start
        int end
        int id
        struct_interval * next

    ctypedef struct_interval interval
    
    cdef struct struct_clusternode:
        int start
        int end
        struct_interval *interval_head
        struct_interval *interval_tail
    
    ctypedef struct_clusternode clusternode
    
    cdef struct struct_clustertree:
        int max_dist
        int min_intervals

        struct_clusternode *root
    
    ctypedef struct_clustertree clustertree
    
    cdef struct struct_treeitr:
        struct_treeitr *next
        struct_clusternode *node
    
    ctypedef struct_treeitr treeitr
    
    clusternode* clusternode_insert(clustertree *tree, clusternode *node, int start, int end, int id)
    clustertree* create_clustertree(int max_dist, int min_intervals)
    treeitr* clusteritr(clustertree *tree)
    void freeclusteritr(treeitr *itr)
    void free_tree(clustertree *tree)

cdef class ClusterTree:
    cdef clustertree *tree
    cdef int mincols
    cdef int minregions
    
    def __cinit__(self, mincols, minregions):
        self.tree = create_clustertree(mincols, minregions)
        self.mincols = mincols
        self.minregions = minregions
        
    def __dealloc__(self):
        free_tree(self.tree)
    
    def insert(self, s, e, id):
        ''' Insert an interval with start, end, id as parameters'''
        if s > e: raise ValueError("Interval start must be before end")
        self.tree.root = clusternode_insert(self.tree, self.tree.root, s, e, id)
                
    def getregions(self):
        ''' Returns a list clusters in ascending order of starting position.
            Each cluster is a tuple of (start, end, [sorted ids of intervals in cluster])
        
        tree = ClusterTree(0, 0)
        Insert (6, 7, 1), (1, 2, 3), (9, 10, 2), (3, 4, 0), (3, 8, 4)
        tree.getregions() returns [(1, 2, [3]), (3, 8, [0, 1, 4]), (9, 10, [2])]
        '''
        cdef treeitr *itr
        cdef interval *ival
        
        regions = []
        itr = clusteritr(self.tree)
        
        while (itr):
            ids = []
            ival = itr.node.interval_head
            while (ival):
                ids.append(ival.id)
                ival = ival.next

            regions.append( (itr.node.start, itr.node.end, sorted(ids)) )
            itr = itr.next
        freeclusteritr(itr)
        return regions
        
    def getlines(self):
        ''' Similar to getregions except it just returns a list of ids of intervals
            The above example would return [3, 0, 1, 4, 2]
         '''
        cdef treeitr *itr
        cdef interval *ival
        
        lines = []
        itr = clusteritr(self.tree)
        
        while (itr):
            ids = []
            ival = itr.node.interval_head
            while (ival):
                ids.append(ival.id)
                ival = ival.next
            
            lines.extend(sorted(ids))
            itr = itr.next
        freeclusteritr(itr)
        return lines
        
