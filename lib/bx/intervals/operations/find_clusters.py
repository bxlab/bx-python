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

import math
import random

from bx.intervals.cluster import ClusterTree
from bx.intervals.io import GenomicInterval


def find_clusters(reader, mincols=1, minregions=2):
    extra = dict()
    chroms = dict()
    linenum = -1
    for interval in reader:
        linenum += 1
        if not isinstance(interval, GenomicInterval):
            extra[linenum] = interval
        else:
            if interval.chrom not in chroms:
                chroms[interval.chrom] = ClusterTree(mincols, minregions)
            try:
                chroms[interval.chrom].insert(interval.start, interval.end, linenum)
            except OverflowError as e:
                try:
                    # This will work only if reader is a NiceReaderWrapper
                    reader.skipped += 1
                    if reader.skipped < 10:
                        reader.skipped_lines.append((reader.linenum, reader.current_line, str(e)))
                except Exception:
                    pass
                continue
    return chroms, extra


# DEPRECATED: Use the ClusterTree in bx.intervals.cluster for this.
# It does the same thing, but is a C implementation.
class ClusterNode:
    def __init__(self, start, end, linenum, mincols, minregions):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = math.ceil((-1.0 / math.log(.5)) * math.log(-1.0 / (random.uniform(0, 1) - 1)))
        self.start = start
        self.end = end
        self.left = None
        self.right = None
        self.lines = [linenum]
        self.mincols = mincols
        self.minregions = minregions

    def insert(self, start, end, linenum):
        if start - self.mincols > self.end:
            # insert to right tree
            if self.right:
                self.right = self.right.insert(start, end, linenum)
            else:
                self.right = ClusterNode(start, end, linenum, self.mincols, self.minregions)
            # rebalance tree
            if self.priority < self.right.priority:
                return self.rotateleft()
        elif end + self.mincols < self.start:
            # insert to left tree
            if self.left:
                self.left = self.left.insert(start, end, linenum)
            else:
                self.left = ClusterNode(start, end, linenum, self.mincols, self.minregions)
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

    def rotateright(self):
        root = self.left
        self.left = self.left.right
        root.right = self
        return root

    def rotateleft(self):
        root = self.right
        self.right = self.right.left
        root.left = self
        return root

    def push_up(self, topnode):
        # Note: this function does not affect heap property
        # Distance method removed for inline, faster?
        distance = max(self.start, topnode.start) - min(self.end, topnode.end)
        if distance <= self.mincols:
            topnode.start = min(self.start, topnode.start)
            topnode.end = max(self.end, topnode.end)
            for linenum in self.lines:
                topnode.lines.append(linenum)
            if self.right:
                return self.right.push_up(topnode)
            if self.left:
                return self.left.push_up(topnode)
            return None
        if self.end < topnode.start and self.right:
            self.right = self.right.push_up(topnode)
        if self.start > topnode.end and self.left:
            self.left = self.left.push_up(topnode)
        return self

    def getintervals(self):
        if self.left:
            yield from self.left.getintervals(self.minregions)
        if len(self.lines) >= self.minregions:
            yield self.start, self.end
        if self.right:
            yield from self.right.getintervals(self.minregions)

    def getlines(self):
        if self.left:
            yield from self.left.getlines()
        if len(self.lines) >= self.minregions:
            yield from self.lines
        if self.right:
            yield from self.right.getlines()
