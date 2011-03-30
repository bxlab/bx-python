"""
Join two sets of intervals using their overlap as the key.  The
intervals MUST be sorted by chrom(lexicographically),
start(arithmetically) and end(arithmetically).  This works by simply
walking through the inputs in O(n) time.
"""

import psyco_full

import math
import traceback
import fileinput
from warnings import warn

from bx.intervals.io import *
from bx.intervals.operations import *
from quicksect import IntervalTree

def join(leftSet, rightSet, mincols=1, leftfill=True, rightfill=True):
    # Read rightSet into memory:
    rightlen = 0
    leftlen = 0
    rightTree = IntervalTree()
    for item in rightSet:
        if isinstance(item, GenomicInterval):
            rightTree.insert( item, rightSet.linenum, item.fields )
            if rightlen == 0: rightlen = item.nfields

    for interval in leftSet:
        if leftlen == 0 and isinstance(interval, GenomicInterval):
            leftlen = interval.nfields
        if not isinstance(interval, GenomicInterval):
            yield interval
        else:
            result = []
            rightTree.intersect( interval, lambda node: result.append( node ) )
            overlap_not_met = 0
            for item in result:
                if item.start in range(interval.start,interval.end+1) and item.end not in range(interval.start,interval.end+1):
                    overlap = interval.end-item.start
                elif item.end in range(interval.start,interval.end+1) and item.start not in range(interval.start,interval.end+1):
                    overlap = item.end-interval.start
                elif item.start in range(interval.start,interval.end+1) and item.end in range(interval.start,interval.end+1):
                    overlap = item.end-item.start
                else:   #the intersecting item's start and end are outside the interval range
                    overlap = interval.end-interval.start
                if overlap < mincols:
                    overlap_not_met += 1
                    continue
                outfields = list(interval)
                map(outfields.append, item.other)
                setattr( item, "visited", True )
                yield outfields
            if (len(result) == 0 or overlap_not_met == len(result)) and rightfill:
                outfields = list(interval)
                for x in range(rightlen): outfields.append(".")
                yield outfields

    if leftfill:
        def report_unvisited( node, results ):
            if not hasattr(node, "visited"):
                results.append( node )
        results = []
        rightTree.traverse( lambda x: report_unvisited( x, results ) )
        for item in results:
            outfields = list()
            for x in range(leftlen): outfields.append(".")
            map(outfields.append, item.other)
            yield outfields


def interval_cmp(a, b):
    interval1 = a[0]
    interval2 = b[0]
    if not (isinstance(interval1, GenomicInterval) and isinstance(interval2, GenomicInterval)):
        return 0
    # Both are intervals
    if interval1.chrom == interval2.chrom:
        center1 = interval1.start + ((interval1.end - interval1.start) / 2)
        center2 = interval2.start + ((interval2.end - interval2.start) / 2)
        return center1 - center2
    else:
        if interval1.chrom > interval2.chrom:
            return 1
        else:
            return -1

    return 0

def findintersect(interval, sortedlist, mincols):
    # find range of intervals that intersect via a binary search
    # find lower bound
    x = len(sortedlist) / 2
    n = int(math.pow(2,math.ceil(math.log(len(sortedlist),2))))

    not_found = True
    not_done = True
    while not_found and not_done:
        n = n / 2
        if n == 0:
            n = 1
            not_done = False
        if x >= len(sortedlist):
            x -= n
        elif x < 0:
            x += n
        else:
            if findoverlap(sortedlist[x][0], interval) >= mincols:
                not_found = False
            else:
                comp = interval_cmp(sortedlist[x], [interval, 0])
                if comp > 0:
                    x -= n
                else:
                    x += n

    print "\t".join(sortedlist[x][0].fields)
    print "not_found = " + str(not_found)
    if not_found:
        return 0,-1

    lowerbound = x
    middlebound = x
    upperbound = x
    while (lowerbound > -1) and (findoverlap(sortedlist[lowerbound-1][0],interval) >= mincols):
        lowerbound -= 1
    while (upperbound+1 < len(sortedlist)) and (findoverlap(sortedlist[upperbound+1][0],interval) >= mincols):
        upperbound += 1

    return lowerbound, upperbound

def findoverlap(a, b):
    # overlapping
    if a.chrom == b.chrom:
        return min(a.end, b.end) - max(a.start, b.start)
    else:
        return 0
