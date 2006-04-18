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
import pkg_resources
pkg_resources.require( "bx-python" )

import psyco_full

import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *


def find_clusters(reader, mincols=1, minregions=2):
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
        # I don't know how many different randoms, or what type of
        # distribution is prefered for this kind of thing.  I'm not a
        # computer scientist.
        self.priority = random.randint(0,500)
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

    def getall( self, minregions ):
        if self.left:
            for start, end in self.left.getall(minregions):
                yield start, end
        if len(self.lines) >= minregions:
            yield self.start, self.end
        if self.right:
            for start, end in self.right.getall(minregions):
                yield start, end
