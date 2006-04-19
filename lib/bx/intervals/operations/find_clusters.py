#!/usr/bin/env python
"""
Find clusters of intervals within a set of intervals.  A cluster is a
group (of size minregions) of intervals within a specific distance (of
mincols) of each other.

Returns Cluster objects, which have a chrom, start, end, and lines (a
list of linenumbers from the original file).  The original can then be
ran through with the linenumbers to extract clustered regions without
disturbing original order, or the clusters may themselves be written
as intervals.
"""

import random
import math
import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *


def find_clusters(reader, mincols=1):
    chroms = dict()
    linenum = -1
    for interval in reader:
        if not type( interval ) is GenomicInterval: continue
        linenum += 1
        if interval.chrom in chroms:
             chroms[interval.chrom] = chroms[interval.chrom].insert(interval.start, interval.end, linenum)
        else:
             chroms[interval.chrom] = ClusterNode(interval.start, interval.end, linenum, mincols)
    return chroms
   
class ClusterNode( object ):
    def __init__( self, start, end, linenum, mincols ):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = math.ceil( (-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
        self.start = start
        self.end = end
        self.left = None
        self.right = None
        self.lines = [linenum]
        self.mincols = mincols
        
    def insert( self, start, end, linenum ):
        if start - self.mincols > self.end:
            # insert to right tree
            if self.right:
                self.right = self.right.insert( start, end, linenum )
            else:
                self.right = ClusterNode(start, end, linenum, self.mincols)
            # rebalance tree
            if self.priority < self.right.priority:
                return self.rotateleft()
        elif end + self.mincols < self.start:
            # insert to left tree
            if self.left:
                self.left = self.left.insert( start, end, linenum )
            else:
                self.left = ClusterNode(start, end, linenum, self.mincols)
            # rebalance tree
            if self.priority < self.left.priority:
                return self.rotateright()
        else:
            # insert here
            self.start = min(self.start, start)
            self.end = max(self.end, end)
            self.lines.append(linenum)
            # recursive call to push nodes up
            if self.left:
                self.left = self.left.push_up(self)
            if self.right:
                self.right = self.right.push_up(self)
        return self

    def rotateright( self ):
        root = self.left
        self.left = self.left.right
        root.right = self
        return root
        
    def rotateleft( self ):
        root = self.right
        self.right = self.right.left
        root.left = self
        return root
        
    def push_up( self, topnode ):
        # Note: this function does not affect heap property
        # Distance method removed for inline, faster?
        distance = 0
        if self.start > topnode.end:
            distance = self.start - topnode.end
        elif topnode.start > self.end:
            distance = topnode.start - self.end
        if distance < self.mincols:
            topnode.start = min(self.start, topnode.start)
            topnode.end = min(self.end, topnode.end)
            for linenum in self.lines:
                topnode.lines.append(linenum)
            if self.right:
                return self.right.push_up( topnode )
            if self.left:
                return self.left.push_up( topnode )
            return None
        if self.end < topnode.start and self.right:
            self.right = self.right.push_up( topnode )
        if self.start > topnode.end and self.left:
            self.left = self.left.push_up( topnode )
        return self

    def getintervals( self, minregions ):
        if self.left:
            for start, end in self.left.getintervals(minregions):
                yield start, end
        if len(self.lines) >= minregions:
            yield self.start, self.end
        if self.right:
            for start, end in self.right.getintervals(minregions):
                yield start, end

    def getlines( self, minregions ):
        if self.left:
            for line in self.left.getlines(minregions):
                yield line
        if len(self.lines) >= minregions:
            for line in self.lines:
                yield line
        if self.right:
            for line in self.right.getlines(minregions):
                yield line
                
## def main():
##     f1 = fileinput.FileInput("big.bed")
##     g1 = GenomicIntervalReader(f1)
##     returntree = find_clusters(g1, mincols=50)
##     for chrom, value in returntree.items():
##         for start, end in value.getintervals(2):
##            print chrom+"\t"+str(start)+"\t"+str(end)
##         for line in value.getlines(2):
##             print "Line:\t"+str(line)

## main()
